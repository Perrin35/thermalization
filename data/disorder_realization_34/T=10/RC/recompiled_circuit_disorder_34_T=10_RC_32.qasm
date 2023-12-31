OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.37378398) q[0];
sx q[0];
rz(-2.7019579) q[0];
sx q[0];
rz(-3.0602732) q[0];
rz(-2.4913139) q[1];
sx q[1];
rz(4.4250017) q[1];
sx q[1];
rz(11.783574) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7711174) q[0];
sx q[0];
rz(-1.3210216) q[0];
sx q[0];
rz(-3.0779548) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15687234) q[2];
sx q[2];
rz(-1.9777159) q[2];
sx q[2];
rz(-3.0770609) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4408735) q[1];
sx q[1];
rz(-2.6881725) q[1];
sx q[1];
rz(-1.2749626) q[1];
rz(-1.4814032) q[3];
sx q[3];
rz(-1.2561744) q[3];
sx q[3];
rz(-0.50953007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.89447442) q[2];
sx q[2];
rz(-1.0049745) q[2];
sx q[2];
rz(-3.0257814) q[2];
rz(1.5995021) q[3];
sx q[3];
rz(-0.096297979) q[3];
sx q[3];
rz(-2.0882864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88749921) q[0];
sx q[0];
rz(-0.54953456) q[0];
sx q[0];
rz(2.9462573) q[0];
rz(0.37503606) q[1];
sx q[1];
rz(-1.476036) q[1];
sx q[1];
rz(-2.9017752) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0175184) q[0];
sx q[0];
rz(-2.8483609) q[0];
sx q[0];
rz(-2.571884) q[0];
rz(1.8717143) q[2];
sx q[2];
rz(-1.7728724) q[2];
sx q[2];
rz(2.2965455) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1344192) q[1];
sx q[1];
rz(-1.6234142) q[1];
sx q[1];
rz(-2.3514071) q[1];
rz(-pi) q[2];
rz(-1.5870985) q[3];
sx q[3];
rz(-1.5566412) q[3];
sx q[3];
rz(-2.6027038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5043162) q[2];
sx q[2];
rz(-0.80233032) q[2];
sx q[2];
rz(1.8117388) q[2];
rz(1.3416393) q[3];
sx q[3];
rz(-1.642671) q[3];
sx q[3];
rz(-1.3247066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.864569) q[0];
sx q[0];
rz(-1.8947911) q[0];
sx q[0];
rz(2.4734316) q[0];
rz(1.6502624) q[1];
sx q[1];
rz(-2.4490093) q[1];
sx q[1];
rz(2.0756762) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50915584) q[0];
sx q[0];
rz(-1.8189438) q[0];
sx q[0];
rz(-3.0993673) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5553603) q[2];
sx q[2];
rz(-1.6931603) q[2];
sx q[2];
rz(0.33081474) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8733858) q[1];
sx q[1];
rz(-1.2899439) q[1];
sx q[1];
rz(-0.98209776) q[1];
x q[2];
rz(3.021574) q[3];
sx q[3];
rz(-0.91306049) q[3];
sx q[3];
rz(2.0002055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9540017) q[2];
sx q[2];
rz(-1.9090586) q[2];
sx q[2];
rz(0.032547396) q[2];
rz(0.35999808) q[3];
sx q[3];
rz(-2.0149752) q[3];
sx q[3];
rz(2.3500197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86768326) q[0];
sx q[0];
rz(-1.5861347) q[0];
sx q[0];
rz(2.4348863) q[0];
rz(1.2061521) q[1];
sx q[1];
rz(-0.34148347) q[1];
sx q[1];
rz(-1.4867841) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.810881) q[0];
sx q[0];
rz(-1.2905621) q[0];
sx q[0];
rz(-2.7964554) q[0];
x q[1];
rz(1.2232259) q[2];
sx q[2];
rz(-0.98993694) q[2];
sx q[2];
rz(-2.6207404) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.023991) q[1];
sx q[1];
rz(-1.7204493) q[1];
sx q[1];
rz(-0.28848044) q[1];
x q[2];
rz(1.5634469) q[3];
sx q[3];
rz(-2.3971933) q[3];
sx q[3];
rz(3.0878382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5965745) q[2];
sx q[2];
rz(-2.2480965) q[2];
sx q[2];
rz(-1.0085227) q[2];
rz(-1.0962983) q[3];
sx q[3];
rz(-1.9129646) q[3];
sx q[3];
rz(0.1299468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1059882) q[0];
sx q[0];
rz(-0.85177079) q[0];
sx q[0];
rz(1.460357) q[0];
rz(1.5885072) q[1];
sx q[1];
rz(-1.2358783) q[1];
sx q[1];
rz(3.1255186) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2747297) q[0];
sx q[0];
rz(-0.99039536) q[0];
sx q[0];
rz(-0.6443364) q[0];
x q[1];
rz(-1.5889421) q[2];
sx q[2];
rz(-2.7108253) q[2];
sx q[2];
rz(-0.22242966) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.88735089) q[1];
sx q[1];
rz(-0.72853959) q[1];
sx q[1];
rz(-0.90708797) q[1];
x q[2];
rz(-3.1297242) q[3];
sx q[3];
rz(-2.7792645) q[3];
sx q[3];
rz(-0.12479898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0323223) q[2];
sx q[2];
rz(-1.0291928) q[2];
sx q[2];
rz(-2.5197022) q[2];
rz(1.0970998) q[3];
sx q[3];
rz(-2.3648839) q[3];
sx q[3];
rz(0.2203075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19875232) q[0];
sx q[0];
rz(-0.0033012882) q[0];
sx q[0];
rz(-0.90676701) q[0];
rz(0.81470195) q[1];
sx q[1];
rz(-2.4532313) q[1];
sx q[1];
rz(-1.2247359) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3186) q[0];
sx q[0];
rz(-2.4791414) q[0];
sx q[0];
rz(-1.9255161) q[0];
rz(-pi) q[1];
rz(-1.9903509) q[2];
sx q[2];
rz(-0.47861368) q[2];
sx q[2];
rz(1.4065557) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8606557) q[1];
sx q[1];
rz(-0.60815629) q[1];
sx q[1];
rz(1.4742736) q[1];
rz(2.8858658) q[3];
sx q[3];
rz(-2.7284107) q[3];
sx q[3];
rz(-2.7030088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1427052) q[2];
sx q[2];
rz(-2.3110516) q[2];
sx q[2];
rz(0.93969807) q[2];
rz(2.9283004) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(-1.7512158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1619103) q[0];
sx q[0];
rz(-0.96452159) q[0];
sx q[0];
rz(-0.58037037) q[0];
rz(-2.0866108) q[1];
sx q[1];
rz(-1.6886728) q[1];
sx q[1];
rz(0.70077983) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38802559) q[0];
sx q[0];
rz(-1.1209079) q[0];
sx q[0];
rz(3.0168424) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31957303) q[2];
sx q[2];
rz(-1.459889) q[2];
sx q[2];
rz(0.061352913) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.80220561) q[1];
sx q[1];
rz(-1.4805668) q[1];
sx q[1];
rz(-1.0046563) q[1];
rz(-pi) q[2];
rz(-1.4768836) q[3];
sx q[3];
rz(-2.4086773) q[3];
sx q[3];
rz(0.35009271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7769527) q[2];
sx q[2];
rz(-0.31034714) q[2];
sx q[2];
rz(3.11943) q[2];
rz(0.74470216) q[3];
sx q[3];
rz(-1.1170758) q[3];
sx q[3];
rz(-2.7409592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78684029) q[0];
sx q[0];
rz(-2.1180034) q[0];
sx q[0];
rz(-1.4165075) q[0];
rz(1.3757061) q[1];
sx q[1];
rz(-1.7388652) q[1];
sx q[1];
rz(0.8917121) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.007908) q[0];
sx q[0];
rz(-2.3500751) q[0];
sx q[0];
rz(-1.3484811) q[0];
rz(-pi) q[1];
rz(0.58015577) q[2];
sx q[2];
rz(-1.9556502) q[2];
sx q[2];
rz(-0.58154026) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5276287) q[1];
sx q[1];
rz(-0.91273897) q[1];
sx q[1];
rz(-1.0440473) q[1];
rz(0.31657747) q[3];
sx q[3];
rz(-0.93571767) q[3];
sx q[3];
rz(-0.77565912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6531758) q[2];
sx q[2];
rz(-1.9501053) q[2];
sx q[2];
rz(-2.3573504) q[2];
rz(-0.50576058) q[3];
sx q[3];
rz(-2.2890746) q[3];
sx q[3];
rz(-2.908356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6655675) q[0];
sx q[0];
rz(-1.5371756) q[0];
sx q[0];
rz(2.419557) q[0];
rz(0.33323914) q[1];
sx q[1];
rz(-1.1958586) q[1];
sx q[1];
rz(-1.3649712) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38395912) q[0];
sx q[0];
rz(-0.74476349) q[0];
sx q[0];
rz(3.0804759) q[0];
rz(-pi) q[1];
rz(1.7609673) q[2];
sx q[2];
rz(-2.724078) q[2];
sx q[2];
rz(-1.6585569) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.078799876) q[1];
sx q[1];
rz(-1.2275057) q[1];
sx q[1];
rz(-0.72453665) q[1];
rz(0.19946675) q[3];
sx q[3];
rz(-1.8623127) q[3];
sx q[3];
rz(0.22649543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1086796) q[2];
sx q[2];
rz(-1.379517) q[2];
sx q[2];
rz(1.9101248) q[2];
rz(3.1075297) q[3];
sx q[3];
rz(-1.27682) q[3];
sx q[3];
rz(-0.6347707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.05474) q[0];
sx q[0];
rz(-0.56607902) q[0];
sx q[0];
rz(-1.4779133) q[0];
rz(1.0832896) q[1];
sx q[1];
rz(-1.3996841) q[1];
sx q[1];
rz(-0.96819425) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6598845) q[0];
sx q[0];
rz(-1.1904926) q[0];
sx q[0];
rz(2.7754521) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2553026) q[2];
sx q[2];
rz(-2.5286525) q[2];
sx q[2];
rz(-2.431589) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.72659513) q[1];
sx q[1];
rz(-1.2508878) q[1];
sx q[1];
rz(-1.725561) q[1];
x q[2];
rz(-1.8929385) q[3];
sx q[3];
rz(-1.570574) q[3];
sx q[3];
rz(-0.9957046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0633462) q[2];
sx q[2];
rz(-2.4137256) q[2];
sx q[2];
rz(-0.001312288) q[2];
rz(-2.0007658) q[3];
sx q[3];
rz(-1.3169378) q[3];
sx q[3];
rz(1.7989981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.39682) q[0];
sx q[0];
rz(-2.0347432) q[0];
sx q[0];
rz(1.9532935) q[0];
rz(2.7753579) q[1];
sx q[1];
rz(-1.2013422) q[1];
sx q[1];
rz(1.3399301) q[1];
rz(2.3920849) q[2];
sx q[2];
rz(-1.3386249) q[2];
sx q[2];
rz(-2.7809536) q[2];
rz(2.7868174) q[3];
sx q[3];
rz(-1.0412024) q[3];
sx q[3];
rz(2.1333221) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
