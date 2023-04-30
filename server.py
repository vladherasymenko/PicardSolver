from flask import Flask, render_template, request
from main import picard_general, list_to_formula


GLOBAL_PRECISION = 3

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    return render_template('input.html')

@app.route('/solve', methods=['GET', 'POST'])
def solve_equations():
    if request.method == 'POST':
        try :
            num_eqns = int(request.form['num_eqns'])
            x0 = float(request.form['x0'])
            func_strings = [request.form[f'eqn{i+1}'] for i in range(num_eqns)]
            init_conds = {}
            for i in range(num_eqns):
                eqn_name = f'y{i+1}'
                eqn_val = float(request.form[f'{eqn_name}_0'])
                init_conds[eqn_name] = [eqn_val]
            taylor_order = 5 #int(request.form['taylor_order'])
            num_iters = int(request.form['iters'])
            solutions = picard_general(func_strings, x0, init_conds, taylor_order, num_iters)
            for key in solutions.keys():
                solutions[key] = list_to_formula(solutions[key], precision=GLOBAL_PRECISION)
            return render_template('output.html', solution=solutions, num_iters=num_iters)
        except:
            return render_template('Error.html')
    else:
        return render_template('input.html')

if __name__ == '__main__':
    app.run(debug=True)