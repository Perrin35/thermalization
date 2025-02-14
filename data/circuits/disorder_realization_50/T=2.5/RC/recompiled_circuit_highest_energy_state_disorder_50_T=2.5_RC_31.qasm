OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7293575) q[0];
sx q[0];
rz(-2.5492302) q[0];
sx q[0];
rz(1.2047729) q[0];
rz(-2.0545948) q[1];
sx q[1];
rz(-0.49538651) q[1];
sx q[1];
rz(-1.187721) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3760494) q[0];
sx q[0];
rz(-2.5693469) q[0];
sx q[0];
rz(-2.8577198) q[0];
rz(-pi) q[1];
rz(-2.4815791) q[2];
sx q[2];
rz(-2.4610991) q[2];
sx q[2];
rz(1.6689491) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.17537847) q[1];
sx q[1];
rz(-1.7433651) q[1];
sx q[1];
rz(-0.047974384) q[1];
x q[2];
rz(0.13924573) q[3];
sx q[3];
rz(-1.2069824) q[3];
sx q[3];
rz(-2.0219959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.54186934) q[2];
sx q[2];
rz(-2.5399127) q[2];
sx q[2];
rz(-2.0987233) q[2];
rz(-3.0964105) q[3];
sx q[3];
rz(-2.9729645) q[3];
sx q[3];
rz(1.6056304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0586108) q[0];
sx q[0];
rz(-3.0284212) q[0];
sx q[0];
rz(0.97852069) q[0];
rz(1.9784031) q[1];
sx q[1];
rz(-0.99449831) q[1];
sx q[1];
rz(1.3226343) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.103687) q[0];
sx q[0];
rz(-0.9390489) q[0];
sx q[0];
rz(0.61236713) q[0];
rz(1.2873994) q[2];
sx q[2];
rz(-1.3958418) q[2];
sx q[2];
rz(2.8307479) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.757431) q[1];
sx q[1];
rz(-2.4520284) q[1];
sx q[1];
rz(-1.227725) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.421881) q[3];
sx q[3];
rz(-2.2595539) q[3];
sx q[3];
rz(-1.0721579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.80403745) q[2];
sx q[2];
rz(-0.20122169) q[2];
sx q[2];
rz(0.03841722) q[2];
rz(-1.5086959) q[3];
sx q[3];
rz(-1.326606) q[3];
sx q[3];
rz(-1.647515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9134193) q[0];
sx q[0];
rz(-2.4689624) q[0];
sx q[0];
rz(-2.9056554) q[0];
rz(1.1783696) q[1];
sx q[1];
rz(-1.447568) q[1];
sx q[1];
rz(-0.49740121) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0097522) q[0];
sx q[0];
rz(-1.3924283) q[0];
sx q[0];
rz(-1.388489) q[0];
x q[1];
rz(-0.8151594) q[2];
sx q[2];
rz(-0.54537702) q[2];
sx q[2];
rz(-2.8653685) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.062089109) q[1];
sx q[1];
rz(-1.5487284) q[1];
sx q[1];
rz(-2.7479321) q[1];
x q[2];
rz(-0.24876066) q[3];
sx q[3];
rz(-1.3205055) q[3];
sx q[3];
rz(0.32969013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1824789) q[2];
sx q[2];
rz(-1.6751869) q[2];
sx q[2];
rz(-1.9100995) q[2];
rz(1.4500827) q[3];
sx q[3];
rz(-1.3030038) q[3];
sx q[3];
rz(2.2836397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0602144) q[0];
sx q[0];
rz(-2.2980818) q[0];
sx q[0];
rz(-0.1928992) q[0];
rz(-1.3554205) q[1];
sx q[1];
rz(-1.5813634) q[1];
sx q[1];
rz(0.17098175) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6714013) q[0];
sx q[0];
rz(-0.51958749) q[0];
sx q[0];
rz(1.6891102) q[0];
x q[1];
rz(2.5582983) q[2];
sx q[2];
rz(-2.3153044) q[2];
sx q[2];
rz(0.14329958) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1384004) q[1];
sx q[1];
rz(-2.698771) q[1];
sx q[1];
rz(-1.3243616) q[1];
x q[2];
rz(-0.90733068) q[3];
sx q[3];
rz(-1.7620834) q[3];
sx q[3];
rz(0.80463791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.96044797) q[2];
sx q[2];
rz(-1.9713216) q[2];
sx q[2];
rz(-2.81874) q[2];
rz(-1.3747831) q[3];
sx q[3];
rz(-2.1026473) q[3];
sx q[3];
rz(0.84097451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5906931) q[0];
sx q[0];
rz(-2.1192079) q[0];
sx q[0];
rz(0.92556959) q[0];
rz(-3.0135221) q[1];
sx q[1];
rz(-1.6431199) q[1];
sx q[1];
rz(-0.099460348) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8680222) q[0];
sx q[0];
rz(-1.5708459) q[0];
sx q[0];
rz(2.2873061) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36014135) q[2];
sx q[2];
rz(-1.8966873) q[2];
sx q[2];
rz(1.2831836) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4872553) q[1];
sx q[1];
rz(-1.0434101) q[1];
sx q[1];
rz(-0.22589639) q[1];
rz(-pi) q[2];
rz(-1.8811536) q[3];
sx q[3];
rz(-2.884293) q[3];
sx q[3];
rz(2.8043487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1818992) q[2];
sx q[2];
rz(-1.5965896) q[2];
sx q[2];
rz(0.37364513) q[2];
rz(-2.2419808) q[3];
sx q[3];
rz(-0.74266946) q[3];
sx q[3];
rz(2.0067748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9099092) q[0];
sx q[0];
rz(-2.9250356) q[0];
sx q[0];
rz(-0.59189558) q[0];
rz(0.20467155) q[1];
sx q[1];
rz(-0.99212956) q[1];
sx q[1];
rz(0.051102292) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6986324) q[0];
sx q[0];
rz(-1.3784096) q[0];
sx q[0];
rz(-1.7257742) q[0];
rz(-pi) q[1];
rz(-1.899753) q[2];
sx q[2];
rz(-0.68495132) q[2];
sx q[2];
rz(-0.16835131) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.40254996) q[1];
sx q[1];
rz(-1.2616757) q[1];
sx q[1];
rz(-1.2039295) q[1];
x q[2];
rz(-2.3971812) q[3];
sx q[3];
rz(-2.4648414) q[3];
sx q[3];
rz(1.6729421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11703141) q[2];
sx q[2];
rz(-2.0270963) q[2];
sx q[2];
rz(1.0943817) q[2];
rz(2.7545605) q[3];
sx q[3];
rz(-0.47457591) q[3];
sx q[3];
rz(-2.835137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7759906) q[0];
sx q[0];
rz(-2.2633573) q[0];
sx q[0];
rz(-2.8940417) q[0];
rz(1.4546825) q[1];
sx q[1];
rz(-1.2431815) q[1];
sx q[1];
rz(-0.036458485) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55900967) q[0];
sx q[0];
rz(-1.8549619) q[0];
sx q[0];
rz(2.3733632) q[0];
rz(-pi) q[1];
rz(0.22561947) q[2];
sx q[2];
rz(-0.7157514) q[2];
sx q[2];
rz(2.4168487) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.79629512) q[1];
sx q[1];
rz(-2.4234001) q[1];
sx q[1];
rz(2.3188306) q[1];
rz(-pi) q[2];
rz(0.074847742) q[3];
sx q[3];
rz(-1.6597431) q[3];
sx q[3];
rz(-1.68881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0869007) q[2];
sx q[2];
rz(-0.96618235) q[2];
sx q[2];
rz(1.3521693) q[2];
rz(-1.2482721) q[3];
sx q[3];
rz(-1.5779481) q[3];
sx q[3];
rz(1.9498391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2317113) q[0];
sx q[0];
rz(-2.8928962) q[0];
sx q[0];
rz(1.4924208) q[0];
rz(-0.17852783) q[1];
sx q[1];
rz(-1.8793722) q[1];
sx q[1];
rz(-2.5915204) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9880542) q[0];
sx q[0];
rz(-0.037153553) q[0];
sx q[0];
rz(-2.1376993) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3230611) q[2];
sx q[2];
rz(-1.5093409) q[2];
sx q[2];
rz(-2.2814192) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4143205) q[1];
sx q[1];
rz(-1.7309233) q[1];
sx q[1];
rz(0.54267197) q[1];
rz(-0.40591235) q[3];
sx q[3];
rz(-1.1391057) q[3];
sx q[3];
rz(-1.5188876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.34755808) q[2];
sx q[2];
rz(-1.4057691) q[2];
sx q[2];
rz(-2.3883635) q[2];
rz(2.8473162) q[3];
sx q[3];
rz(-0.080065057) q[3];
sx q[3];
rz(-0.068597138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2733961) q[0];
sx q[0];
rz(-0.28565872) q[0];
sx q[0];
rz(-1.4519325) q[0];
rz(-2.946335) q[1];
sx q[1];
rz(-2.0611019) q[1];
sx q[1];
rz(1.6498227) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0367835) q[0];
sx q[0];
rz(-2.0441136) q[0];
sx q[0];
rz(1.9154873) q[0];
rz(-3.0259616) q[2];
sx q[2];
rz(-2.3732269) q[2];
sx q[2];
rz(2.7344784) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50290132) q[1];
sx q[1];
rz(-1.4152621) q[1];
sx q[1];
rz(-3.0763094) q[1];
x q[2];
rz(1.5774861) q[3];
sx q[3];
rz(-1.5433528) q[3];
sx q[3];
rz(-2.5344078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1008272) q[2];
sx q[2];
rz(-2.6746174) q[2];
sx q[2];
rz(0.85774285) q[2];
rz(-1.4486676) q[3];
sx q[3];
rz(-1.492615) q[3];
sx q[3];
rz(-2.9878476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1345074) q[0];
sx q[0];
rz(-2.9869098) q[0];
sx q[0];
rz(-0.78306985) q[0];
rz(-2.018441) q[1];
sx q[1];
rz(-1.4479366) q[1];
sx q[1];
rz(0.060308594) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9911847) q[0];
sx q[0];
rz(-1.6643608) q[0];
sx q[0];
rz(1.9267328) q[0];
rz(-pi) q[1];
rz(-0.98452703) q[2];
sx q[2];
rz(-1.6590377) q[2];
sx q[2];
rz(-0.72438223) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7544514) q[1];
sx q[1];
rz(-1.5031414) q[1];
sx q[1];
rz(0.67004244) q[1];
rz(-pi) q[2];
rz(-2.8496212) q[3];
sx q[3];
rz(-0.37759104) q[3];
sx q[3];
rz(2.1462314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.15410885) q[2];
sx q[2];
rz(-2.2769112) q[2];
sx q[2];
rz(1.8314499) q[2];
rz(1.3335258) q[3];
sx q[3];
rz(-1.3430877) q[3];
sx q[3];
rz(1.8346627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6780728) q[0];
sx q[0];
rz(-2.0250043) q[0];
sx q[0];
rz(2.9153839) q[0];
rz(2.3851867) q[1];
sx q[1];
rz(-0.49497985) q[1];
sx q[1];
rz(2.3351647) q[1];
rz(1.5802154) q[2];
sx q[2];
rz(-1.1404372) q[2];
sx q[2];
rz(-2.5802317) q[2];
rz(-1.2213094) q[3];
sx q[3];
rz(-1.3987473) q[3];
sx q[3];
rz(-1.7940298) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
