OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7678087) q[0];
sx q[0];
rz(5.8435506) q[0];
sx q[0];
rz(6.2018659) q[0];
rz(-2.4913139) q[1];
sx q[1];
rz(-1.8581837) q[1];
sx q[1];
rz(2.3587956) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7711174) q[0];
sx q[0];
rz(-1.3210216) q[0];
sx q[0];
rz(-0.063637861) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2230258) q[2];
sx q[2];
rz(-0.43453056) q[2];
sx q[2];
rz(0.44473106) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.767747) q[1];
sx q[1];
rz(-1.1384374) q[1];
sx q[1];
rz(0.14111622) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26773914) q[3];
sx q[3];
rz(-0.32666884) q[3];
sx q[3];
rz(-0.79145811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2471182) q[2];
sx q[2];
rz(-1.0049745) q[2];
sx q[2];
rz(3.0257814) q[2];
rz(-1.5420906) q[3];
sx q[3];
rz(-0.096297979) q[3];
sx q[3];
rz(1.0533062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88749921) q[0];
sx q[0];
rz(-2.5920581) q[0];
sx q[0];
rz(-0.19533531) q[0];
rz(-0.37503606) q[1];
sx q[1];
rz(-1.476036) q[1];
sx q[1];
rz(2.9017752) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.996802) q[0];
sx q[0];
rz(-1.7273434) q[0];
sx q[0];
rz(-0.24897225) q[0];
rz(1.8717143) q[2];
sx q[2];
rz(-1.7728724) q[2];
sx q[2];
rz(-0.84504715) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6310196) q[1];
sx q[1];
rz(-0.78200713) q[1];
sx q[1];
rz(1.4960947) q[1];
rz(-0.85571839) q[3];
sx q[3];
rz(-3.120003) q[3];
sx q[3];
rz(1.3947226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5043162) q[2];
sx q[2];
rz(-2.3392623) q[2];
sx q[2];
rz(-1.8117388) q[2];
rz(-1.7999533) q[3];
sx q[3];
rz(-1.4989217) q[3];
sx q[3];
rz(1.3247066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.864569) q[0];
sx q[0];
rz(-1.2468015) q[0];
sx q[0];
rz(-2.4734316) q[0];
rz(1.6502624) q[1];
sx q[1];
rz(-2.4490093) q[1];
sx q[1];
rz(2.0756762) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50915584) q[0];
sx q[0];
rz(-1.3226489) q[0];
sx q[0];
rz(-3.0993673) q[0];
rz(-pi) q[1];
rz(-0.22914825) q[2];
sx q[2];
rz(-2.5742968) q[2];
sx q[2];
rz(2.0958401) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.26820688) q[1];
sx q[1];
rz(-1.2899439) q[1];
sx q[1];
rz(0.98209776) q[1];
rz(-pi) q[2];
rz(1.4170253) q[3];
sx q[3];
rz(-0.66699281) q[3];
sx q[3];
rz(-2.1949777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.187591) q[2];
sx q[2];
rz(-1.9090586) q[2];
sx q[2];
rz(-3.1090453) q[2];
rz(-0.35999808) q[3];
sx q[3];
rz(-2.0149752) q[3];
sx q[3];
rz(0.79157296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2739094) q[0];
sx q[0];
rz(-1.5554579) q[0];
sx q[0];
rz(-2.4348863) q[0];
rz(-1.9354405) q[1];
sx q[1];
rz(-0.34148347) q[1];
sx q[1];
rz(1.6548086) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1409722) q[0];
sx q[0];
rz(-1.9019433) q[0];
sx q[0];
rz(1.2739869) q[0];
rz(-pi) q[1];
rz(-1.2232259) q[2];
sx q[2];
rz(-2.1516557) q[2];
sx q[2];
rz(0.52085224) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0123803) q[1];
sx q[1];
rz(-0.32402363) q[1];
sx q[1];
rz(-0.48735168) q[1];
rz(-pi) q[2];
rz(1.5781457) q[3];
sx q[3];
rz(-0.74439936) q[3];
sx q[3];
rz(3.0878382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5450181) q[2];
sx q[2];
rz(-0.89349616) q[2];
sx q[2];
rz(-1.0085227) q[2];
rz(2.0452943) q[3];
sx q[3];
rz(-1.228628) q[3];
sx q[3];
rz(-0.1299468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1059882) q[0];
sx q[0];
rz(-0.85177079) q[0];
sx q[0];
rz(1.6812356) q[0];
rz(1.5885072) q[1];
sx q[1];
rz(-1.9057143) q[1];
sx q[1];
rz(-3.1255186) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0468633) q[0];
sx q[0];
rz(-1.0445347) q[0];
sx q[0];
rz(2.2577283) q[0];
rz(-pi) q[1];
rz(0.0083382567) q[2];
sx q[2];
rz(-1.1401046) q[2];
sx q[2];
rz(-2.8991933) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2542418) q[1];
sx q[1];
rz(-2.4130531) q[1];
sx q[1];
rz(-0.90708797) q[1];
rz(-0.36230476) q[3];
sx q[3];
rz(-1.575003) q[3];
sx q[3];
rz(1.6844974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1092704) q[2];
sx q[2];
rz(-2.1123999) q[2];
sx q[2];
rz(-2.5197022) q[2];
rz(2.0444929) q[3];
sx q[3];
rz(-0.77670875) q[3];
sx q[3];
rz(-2.9212852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19875232) q[0];
sx q[0];
rz(-3.1382914) q[0];
sx q[0];
rz(2.2348256) q[0];
rz(-0.81470195) q[1];
sx q[1];
rz(-2.4532313) q[1];
sx q[1];
rz(-1.9168568) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46366102) q[0];
sx q[0];
rz(-1.3555129) q[0];
sx q[0];
rz(-2.2023375) q[0];
rz(2.9333026) q[2];
sx q[2];
rz(-2.0049094) q[2];
sx q[2];
rz(-1.2693894) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9310301) q[1];
sx q[1];
rz(-1.6258874) q[1];
sx q[1];
rz(-2.1767666) q[1];
rz(0.25572689) q[3];
sx q[3];
rz(-2.7284107) q[3];
sx q[3];
rz(2.7030088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.99888745) q[2];
sx q[2];
rz(-2.3110516) q[2];
sx q[2];
rz(2.2018946) q[2];
rz(-2.9283004) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(1.7512158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9796824) q[0];
sx q[0];
rz(-2.1770711) q[0];
sx q[0];
rz(-0.58037037) q[0];
rz(-2.0866108) q[1];
sx q[1];
rz(-1.4529198) q[1];
sx q[1];
rz(2.4408128) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0133007) q[0];
sx q[0];
rz(-1.4585146) q[0];
sx q[0];
rz(-2.0237472) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31957303) q[2];
sx q[2];
rz(-1.459889) q[2];
sx q[2];
rz(3.0802397) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82579457) q[1];
sx q[1];
rz(-1.0072395) q[1];
sx q[1];
rz(3.0347996) q[1];
rz(-2.3015162) q[3];
sx q[3];
rz(-1.6335765) q[3];
sx q[3];
rz(-1.9907794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7769527) q[2];
sx q[2];
rz(-2.8312455) q[2];
sx q[2];
rz(-3.11943) q[2];
rz(-0.74470216) q[3];
sx q[3];
rz(-2.0245168) q[3];
sx q[3];
rz(-2.7409592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78684029) q[0];
sx q[0];
rz(-1.0235893) q[0];
sx q[0];
rz(1.4165075) q[0];
rz(1.3757061) q[1];
sx q[1];
rz(-1.4027275) q[1];
sx q[1];
rz(2.2498806) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4211593) q[0];
sx q[0];
rz(-1.4132858) q[0];
sx q[0];
rz(-2.349855) q[0];
x q[1];
rz(-0.58015577) q[2];
sx q[2];
rz(-1.1859425) q[2];
sx q[2];
rz(-0.58154026) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5276287) q[1];
sx q[1];
rz(-0.91273897) q[1];
sx q[1];
rz(-1.0440473) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31657747) q[3];
sx q[3];
rz(-0.93571767) q[3];
sx q[3];
rz(0.77565912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4884168) q[2];
sx q[2];
rz(-1.1914873) q[2];
sx q[2];
rz(-0.78424224) q[2];
rz(-2.6358321) q[3];
sx q[3];
rz(-0.85251802) q[3];
sx q[3];
rz(0.23323664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4760251) q[0];
sx q[0];
rz(-1.5371756) q[0];
sx q[0];
rz(-2.419557) q[0];
rz(0.33323914) q[1];
sx q[1];
rz(-1.9457341) q[1];
sx q[1];
rz(1.3649712) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38395912) q[0];
sx q[0];
rz(-2.3968292) q[0];
sx q[0];
rz(0.061116771) q[0];
rz(1.981609) q[2];
sx q[2];
rz(-1.4940726) q[2];
sx q[2];
rz(3.0551747) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8553798) q[1];
sx q[1];
rz(-0.78821048) q[1];
sx q[1];
rz(2.6469995) q[1];
rz(-pi) q[2];
rz(2.9421259) q[3];
sx q[3];
rz(-1.8623127) q[3];
sx q[3];
rz(2.9150972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1086796) q[2];
sx q[2];
rz(-1.379517) q[2];
sx q[2];
rz(-1.2314679) q[2];
rz(-3.1075297) q[3];
sx q[3];
rz(-1.8647727) q[3];
sx q[3];
rz(2.506822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0868527) q[0];
sx q[0];
rz(-0.56607902) q[0];
sx q[0];
rz(-1.4779133) q[0];
rz(1.0832896) q[1];
sx q[1];
rz(-1.3996841) q[1];
sx q[1];
rz(-0.96819425) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23046872) q[0];
sx q[0];
rz(-1.9096806) q[0];
sx q[0];
rz(1.1662657) q[0];
rz(0.88629006) q[2];
sx q[2];
rz(-2.5286525) q[2];
sx q[2];
rz(-2.431589) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.248368) q[1];
sx q[1];
rz(-1.7176504) q[1];
sx q[1];
rz(-2.8180772) q[1];
x q[2];
rz(1.5700941) q[3];
sx q[3];
rz(-2.8194504) q[3];
sx q[3];
rz(2.567167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0782464) q[2];
sx q[2];
rz(-0.72786704) q[2];
sx q[2];
rz(-0.001312288) q[2];
rz(2.0007658) q[3];
sx q[3];
rz(-1.3169378) q[3];
sx q[3];
rz(1.3425945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7447727) q[0];
sx q[0];
rz(-2.0347432) q[0];
sx q[0];
rz(1.9532935) q[0];
rz(2.7753579) q[1];
sx q[1];
rz(-1.2013422) q[1];
sx q[1];
rz(1.3399301) q[1];
rz(2.8075519) q[2];
sx q[2];
rz(-0.77790778) q[2];
sx q[2];
rz(1.6890656) q[2];
rz(-0.35477521) q[3];
sx q[3];
rz(-1.0412024) q[3];
sx q[3];
rz(2.1333221) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];