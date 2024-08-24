<!DOCTYPE html>
<!---Coding By CodingLab | www.codinglabweb.com--->
<html lang="en">
  <?php $testError = " "; ?>
  <head>
    <meta charset="UTF-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0" />
    <meta http-equiv="X-UA-Compatible" content="ie=edge" />
    <!--<title>Registration Form in HTML CSS</title>-->
    <!---Custom CSS File--->
    <link rel="stylesheet" href="styles.css" />
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
  </head>
  <script>
    function onClickOfBut() {
      document.getElementById("LoadingMessage").innerText = "Please, wait a second ...";

    }
  </script>
  <body>
    <section class="container">
      <header>Multigrid Inputs</header>
      <form action="home.php" class="form" method="post" onsubmit="onClickOfBut()">
        <div class="input-box">
          <label>Insert a valid N:</label>
          <input type="number" name="N" placeholder="Enter N value" required />
        </div>

        <div class="input-box">
          <label>Insert variable a of Poisson problem value</label>
          <input type="text" name="a" placeholder="Enter a" required />
        </div>

        <div class="input-box">
          <label>Insert width</label>
          <input type="text" name="width" placeholder="Enter width value" required />
        </div>

        <div class="column">
          <div class="input-box">
            <label>Multigrid deepness Level</label>
            <input type="number" name="level" placeholder="Enter multigrid level" required />
          </div>
          <div class="input-box">
            <label>Test case</label>
            <div class="select-box">
              <select name="test_selection">
                <option hidden>Tests</option>
                <option value="1">f = -5.0 * exp(x) * exp(-2.0 * y) ;
                        g = exp(x) * exp(-2.0 * y)
                </option>
                <option value="3">f = 1.0 ; g = 0.0</option>
              </select>
            </div>
          </div>
        </div>
        <div class="gender-box">
          <h3>Smoother:</h3>
          <div class="gender-option">
            <div class="gender">
              <input type="radio" id="check-male" name="smoother" value="1" checked />
              <label for="check-male">Jacobi</label>
            </div>
            <div class="gender">
              <input type="radio" id="check-female" name="smoother" value="0" />
              <label for="check-female">Gauss Siedel</label>
            </div>
            <div class="gender">
              <input type="radio" id="check-other" name="smoother" value="2" />
              <label for="check-other">CG</label>
            </div>
          </div>
        </div>
        <button>Submit</button>
        
      </form>
      <?php
        $output = null;  
        
        if (isset($_POST["N"]))
        {
          $N = (int)$_POST["N"];
          $l = (int)$_POST["level"];
          for($temp = 1; $temp < $l; $temp++ )
          {
              $N = $N * 2 - 1;
          }
          $a = $_POST["a"];
          $w = $_POST["width"];
          $test = $_POST["test_selection"];
          $smoot = $_POST["smoother"];
          $output = shell_exec("./Multigrid -n $N -a $a -w $w -ml $l -test $test -smt $smoot");

          if (stripos($output, 'Error:') !== false)
          {
            
            $testError = substr($output, stripos($output, 'Error:'));;
            $output = null;
          }else{
            $testError = " ";
          }
            
        }
        
        
        if ($output != null)
        {
            echo '<canvas  id="myChart" width="400" height="200"></canvas>';
            echo '<a href="x.mtx" download="x.mtx" class="download-button">
                Download solution file
            </a>';

            echo '<a href="MGGS4.txt" download="MGGS4.txt" class="download-button">
                Download convergence file
            </a>';
        }
        
      ?>
      <h4 id="LoadingMessage" style="color: red;"><?php echo $testError; ?></h4>
    </section>
    
    <script>
        // Function to fetch and parse the data
        async function fetchData() {
            const response = await fetch('MGGS4.txt'); // Fetch the data file
            const text = await response.text(); // Get the file content as text

            // Parse the text into an array of numbers
            const values = text.trim().split('\n').map(line => parseFloat(line));

            return values;
        }

        // Function to create the chart
        async function createChart() {
            const values = await fetchData();

            const labels = values.map((_, index) => index + 1);

            // Get the context of the canvas element
            const ctx = document.getElementById('myChart').getContext('2d');

            // Create a new Chart instance
            new Chart(ctx, {
                type: 'line', // Type of chart
                data: {
                    labels: labels, // X-axis labels
                    datasets: [{
                        label: 'Convergence chart', // Legend label
                        data: values, // Data points
                        borderColor: 'rgba(75, 192, 192, 1)', // Line color
                        borderWidth: 1, // Line width
                        fill: false // Disable filling under the line
                    }]
                },
                options: {
                    scales: {
                        x: {
                            title: {
                                display: true,
                                text: 'Iteration-th'
                            }
                        },
                        y: {
                            beginAtZero: false, // Allow Y-axis to adjust based on data
                            title: {
                                display: true,
                                text: 'Value'
                            },
                            ticks: {
                                callback: function(value) {
                                    // Return the value as a string to avoid rounding
                                    return value.toString();
                                }
                            }
                        }
                    }
                }
            });
        }

        // Create the chart when the page loads
        createChart();
    </script>
  </body>
</html>