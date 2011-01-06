#!/usr/bin/php
<?php

if ($_SERVER['argc'] < 2)
{
	$cmd = $_SERVER['argv'][0];
	echo "No basis set specified\n";
	echo "Usage: {$cmd} basis [format]\n";
	exit(1);
}
$basis_name = $_SERVER['argv'][1];
$format = 'Turbomole';
if ($_SERVER['argc'] > 2)
	$format = $_SERVER['argv'][2];

$elements = array('H', 'He', 'Li');

$known_formats = array(
	'Gaussian94' => 'gaussian94',
	'GAMESS-US' => 'gamessus',
	'GAMESS-UK' => 'gamessuk',
	'Turbomole' => 'turbomole',
	'TX93' => 'tx93',
	'Molpro' => 'molpro',
	'MolproInt' => 'molproint',
	'Hondo' => 'hondo',
	'SuperMolecule' => 'supermolecule',
	'Molcas' => 'molcas',
	'HyperChem' => 'hyperchem',
	'Dalton' => 'dalton',
	'deMon-KS' => 'demonks',
	'deMon2k' => 'demon2k',
	'AcesII' => 'aces2'
);
if (!array_key_exists($format, $known_formats))
{
	echo "Unknown output format {$format}\n";
	exit(1);
}

$jar = 'cookies.dat';
$url = 'https://bse.pnl.gov:443/bse/portal/user/anon/js_peid/11535052407933/panel/Main/template/content';

$curl = curl_init($url);
curl_setopt($curl, CURLOPT_COOKIEJAR, $jar);
curl_setopt($curl, CURLOPT_RETURNTRANSFER, true);
$data = curl_exec($curl);

unset($matches);
if (!preg_match_all('/basisSets\[\d+\] = new basisSet\("(.*?)"\)/', $data,
	$matches))
{
	echo "Basis sets not found\n";
	exit(1);
}
$sets = array();
foreach ($matches[1] as $match)
{
	$fields = preg_split('/", *"/', $match);
	$sets[] = $fields;
}

$sel = array();
foreach ($sets as $i => $set)
{
	if ($set[1] == $basis_name)
	        $sel[] = $i;
}
if (!$sel)
{
	echo "No basis sets with name {$basis_name} found\n";
	exit(1);
}
if (count($sel) > 1)
{
	echo "Multiple basis sets with name {$basis_name} found\n";
	exit(1);
}
$set = $sets[reset($sel)];

$url = 'https://bse.pnl.gov:443/bse/portal/user/anon/js_peid/11535052407933/action/portlets.BasisSetAction/template/courier_content/panel/Main/eventSubmit_doDownload/true';
$args = array(
	'bsurl' => $set[0],
	'bsname' => $basis_name,
	'elts' => implode(' ', $elements),
	'format' => $format,
	'minimize' => ''
);
$q = http_build_query($args);
$url .= '?' . $q;

$curl = curl_init($url);
curl_setopt($curl, CURLOPT_COOKIEJAR, $jar);
curl_setopt($curl, CURLOPT_RETURNTRANSFER, true);
//curl_setopt($curl, CURLINFO_HEADER_OUT, true);
curl_setopt($curl, CURLOPT_FOLLOWLOCATION, true);
$data = curl_exec($curl);

unset($matches);
if (preg_match('!<pre .*?>(.*?)</pre>!s', $data, $matches))
{
	$ext = $known_formats[$format];
	file_put_contents("{$basis_name}.{$ext}", $matches[1]);
}
else
{
	echo "Failed to find basis set data\n";
	exit(1);
}

?>
