from pprint import pprint

from pyomo.core import ConcreteModel, Var, Binary, Objective, minimize, ConstraintList
from pyomo.environ import SolverFactory, SolverManagerFactory

days = ['Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat', 'Sun']

shifts = ['morning', 'evening', 'night']
days_shifts = {day: shifts for day in days}

workers = [f'W{i}' for i in range(1, 11)]

model = ConcreteModel()

# Now we can add binary variables representing states of the model.
# Note that we use tuples that will define indexing of the variable, so that later we can access its state for a given
# index using dictionary type access some_variable[index1, index2, index3.....]

# Binary variables representing if given worker on day at shift is scheduled to work.
# We use (worker, day, shift) = 0/1 tuples to represent this information.
model.works = Var(
    ((worker, day, shift) for worker in workers for day in days for shift in days_shifts[day]),
    within=Binary,
    initialize=0
)

# Binary variables representing if worker is necessary - should we use this worker in the workforce.
model.needed = Var(workers, within=Binary, initialize=0)

# Binary variables representing if worker worked on Sunday but not on Saturday - we want to avoid that
# 1 if worker does not work on Sunday but works on Saturday.
model.no_pref = Var(workers, within=Binary, initialize=0)

# Objective function:
def obj_rule(m):
    c = len(workers)
    # Our primary objective is to minimize number of workers needed -> sum(c * m.needed[worker])
    # We multiply by a constant 'c' to make sure that this part of the objective is most important.
    # Secondary objective is no_pref (avoiding scheduling 1 days in the weekends)
    return sum(m.no_pref[worker] for worker in workers) + sum(c * m.needed[worker] for worker in workers)

# Now we can add the objective function to the model and set it to be minimized.
# The objective now is to find a schedule minimizing the number of workers needed and once that is done, also reduce
# the number of workers who have to work on Sundays but not on Saturdays.
# The constant used to multiply needed workers makes sure that this is the primary objective,
model.obj = Objective(rule=obj_rule, sense=minimize)

# Now we can add the constraints that describe our food store.

# We create a set of constraints on the model.
model.constraints = ConstraintList()

# 1. Constraint to make sure that all shifts are assigned and appropriate number of workers are working,
for day in days:
    for shift in days_shifts[day]:
        if day != 'Sunday' and shift in ['morning', 'evening']:
            # Weekdays and Saturday, morning and evenings have 2 workers.
            # Note that constraints are booleans!
            # Here we make sure that for each variable in works for a given worker, day, shift the sum of Binary values
            # is exactly 2.
            model.constraints.add(
                sum(model.works[worker, day, shift] for worker in workers) == 2
            )
        else:
            # For Sundays or night shifts we need only 1 worker.
            model.constraints.add(
                sum(model.works[worker, day, shift] for worker in workers) == 1
            )

# 2. Constraint to make sure each worker does maximally 40 hours of work.
# We are multiplying number of binaries representing that worker works in a given shift by 8 to calculate total hours
# worked.
for worker in workers:
    model.constraints.add(
        sum(model.works[worker, day, shift] * 8 for day in days for shift in days_shifts[day]) <= 40
    )

# 3. Rest between two shifts is greater or equal to 12 hours. Since shift is 8hrs, at least 2 shifts rest is needed
# meaning that only one shift can be worked in any 3.
for worker in workers:
    for day_number in range(len(days)):
        # If starting in the morning, we cannot work any other shift on that day
        model.constraints.add(
            sum(model.works[worker, days[day_number], shift] for shift in days_shifts[days[day_number]]) <= 1
        )

        # If working in the evening, we go to the next day and make sure that we do not work in the morning or evening.
        # Also we have a problem that after Sunday comes Monday - so we have to be able to cycle. We will use modulo
        # operator for that.
        model.constraints.add(
            sum(model.works[worker, days[day_number], shift] for shift in ['evening', 'night']) +
            model.works[worker, days[(day_number + 1) % 7], 'morning'] <= 1
        )

        # If working in the evening we make sure we do not work until next evening. Again the day changes so we have to
        # use modulo operator to cycle.
        model.constraints.add(
            model.works[worker, days[day_number], 'night'] +
            sum(model.works[worker, days[(day_number + 1) % 7], shift] for shift in ['morning', 'evening']) <= 1
        )

# 4. Constraint defining needed number of workers for the objective function. If any worker does any shift they are
# needed.
for worker in workers:
    model.constraints.add(
        sum(model.works[worker, day, shift] for day in days for shift in days_shifts[day]) <=
        10000 * model.needed[worker]
    )

# 5. Constraint defining preference for weekend work.
for worker in workers:
    model.constraints.add(
        sum(model.works[worker, 'Sat', shift] for shift in days_shifts['Sat']) -
        sum(model.works[worker, 'Sun', shift] for shift in days_shifts['Sun']) <= model.no_pref[worker]

    )

# We will use coin or branch and cut solver
opt = SolverFactory('cbc')

# We will use neos server to solve
solver_manager = SolverManagerFactory('neos')

results = solver_manager.solve(model, opt=opt)

def get_workers_needed(needed):
    """Extract to a list the needed workers for the optimal solution."""
    workers_needed = []
    for worker in workers:
        if needed[worker].value == 1:
            workers_needed.append(worker)
    return workers_needed


def get_work_table(works):
    """Build a timetable of the week as a dictionary from the model's optimal solution."""
    week_table = {day: {shift: [] for shift in days_shifts[day]} for day in days}
    for worker in workers:
        for day in days:
            for shift in days_shifts[day]:
                    if works[worker, day, shift].value == 1:
                        week_table[day][shift].append(worker)
    return week_table


def get_no_preference(no_pref):
    """Extract to a list the workers not satisfied with their weekend preference."""
    return [worker for worker in workers if no_pref[worker].value == 1]


workers_needed = get_workers_needed(model.needed)  # dict with the optimal timetable
week_table = get_work_table(model.works)  # list with the required workers
workers_no_pref = get_no_preference(model.no_pref)  # list with the non-satisfied workers (work on Sat but not on Sun)

pprint(workers_needed)
pprint(week_table)
pprint(workers_no_pref)
