<?php

// File path
$filePath = '../src/utilities.cpp';

// Initialize variables
$startReading = false;
$lines = [];

try {
    // Open the file
    $file = new SplFileObject($filePath);

    // Read each line
    while (!$file->eof()) {
        $line = $file->fgets();

        // Check for the start marker
        if (strpos($line, '// FFF') !== false) {
            $startReading = true;
            continue;  // Skip the line containing // FFF
        }

        // Check for the end marker
        if (strpos($line, '// END') !== false) {
            break;  // Stop reading when we reach // END
        }

        // Collect lines between the markers
        if ($startReading) {
            $lines[] = $line;
        }
    }

    // Optionally close the file (handled automatically by PHP when using SplFileObject)
    $file = null;

    // Print or process the captured lines

    $allFunctions = [];
    $count = 0;
    foreach ($lines as $capturedLine) {
        $capturedLine = substr($capturedLine, stripos($capturedLine, 'return') + 6);
        $capturedLine = substr($capturedLine, 0, -(strlen($capturedLine) - stripos($capturedLine, ';')));
        if($count %2 == 0){
            array_push($allFunctions, 'value="' . $count /2 . '">f = ' . $capturedLine . " ; ");
        }else{
            array_push($allFunctions, "g = " . $capturedLine);
        }
        $count++;
    }

} catch (Exception $e) {
    echo 'Error to find functions: ',  $e->getMessage(), "\n";
}

?>
