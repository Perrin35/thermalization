OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3671626) q[0];
sx q[0];
rz(4.055152) q[0];
sx q[0];
rz(11.154296) q[0];
rz(-2.9867759) q[1];
sx q[1];
rz(5.6875416) q[1];
sx q[1];
rz(7.7653801) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1971561) q[0];
sx q[0];
rz(-1.7261788) q[0];
sx q[0];
rz(2.7128501) q[0];
rz(-2.6719195) q[2];
sx q[2];
rz(-2.854752) q[2];
sx q[2];
rz(-1.4925721) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2130148) q[1];
sx q[1];
rz(-1.1482114) q[1];
sx q[1];
rz(-1.0282474) q[1];
x q[2];
rz(1.0079908) q[3];
sx q[3];
rz(-1.4853012) q[3];
sx q[3];
rz(2.8443955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.98510629) q[2];
sx q[2];
rz(-2.6323695) q[2];
sx q[2];
rz(-0.86581725) q[2];
rz(2.1872897) q[3];
sx q[3];
rz(-1.538397) q[3];
sx q[3];
rz(-1.8538063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99825478) q[0];
sx q[0];
rz(-1.7049494) q[0];
sx q[0];
rz(0.026219333) q[0];
rz(1.5401309) q[1];
sx q[1];
rz(-1.5427579) q[1];
sx q[1];
rz(0.96347934) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026982633) q[0];
sx q[0];
rz(-0.61404213) q[0];
sx q[0];
rz(1.5747889) q[0];
rz(-pi) q[1];
rz(-1.0785525) q[2];
sx q[2];
rz(-0.62765593) q[2];
sx q[2];
rz(-0.17222675) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0986833) q[1];
sx q[1];
rz(-0.60815647) q[1];
sx q[1];
rz(-2.152918) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96879949) q[3];
sx q[3];
rz(-2.7235944) q[3];
sx q[3];
rz(2.9434162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5144689) q[2];
sx q[2];
rz(-2.0141979) q[2];
sx q[2];
rz(-0.13452402) q[2];
rz(-2.3965805) q[3];
sx q[3];
rz(-2.9146505) q[3];
sx q[3];
rz(-0.9427332) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9298252) q[0];
sx q[0];
rz(-2.7524502) q[0];
sx q[0];
rz(-2.3441558) q[0];
rz(2.0939317) q[1];
sx q[1];
rz(-2.9918549) q[1];
sx q[1];
rz(-2.581596) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3968351) q[0];
sx q[0];
rz(-2.7042537) q[0];
sx q[0];
rz(1.0768946) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0736018) q[2];
sx q[2];
rz(-2.8376841) q[2];
sx q[2];
rz(1.7169203) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9251688) q[1];
sx q[1];
rz(-0.41453002) q[1];
sx q[1];
rz(0.51149909) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6651354) q[3];
sx q[3];
rz(-1.1343079) q[3];
sx q[3];
rz(2.1863294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3893163) q[2];
sx q[2];
rz(-1.9439149) q[2];
sx q[2];
rz(-2.9690572) q[2];
rz(0.98207384) q[3];
sx q[3];
rz(-1.7445824) q[3];
sx q[3];
rz(-1.0579695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.003222) q[0];
sx q[0];
rz(-1.0976185) q[0];
sx q[0];
rz(0.28451434) q[0];
rz(-2.8248887) q[1];
sx q[1];
rz(-2.7088294) q[1];
sx q[1];
rz(1.8428615) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38628681) q[0];
sx q[0];
rz(-2.5884429) q[0];
sx q[0];
rz(-2.0095216) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1403055) q[2];
sx q[2];
rz(-0.80438559) q[2];
sx q[2];
rz(3.121701) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.87151566) q[1];
sx q[1];
rz(-1.75711) q[1];
sx q[1];
rz(-0.99888505) q[1];
rz(-1.2508568) q[3];
sx q[3];
rz(-2.7971929) q[3];
sx q[3];
rz(-2.875945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5359042) q[2];
sx q[2];
rz(-0.19583344) q[2];
sx q[2];
rz(2.7569125) q[2];
rz(-2.3875333) q[3];
sx q[3];
rz(-1.056517) q[3];
sx q[3];
rz(1.6872905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5383179) q[0];
sx q[0];
rz(-2.1544927) q[0];
sx q[0];
rz(1.3866562) q[0];
rz(-0.23100135) q[1];
sx q[1];
rz(-1.341154) q[1];
sx q[1];
rz(2.8447661) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7291527) q[0];
sx q[0];
rz(-1.531633) q[0];
sx q[0];
rz(-0.57106437) q[0];
x q[1];
rz(-2.3987531) q[2];
sx q[2];
rz(-1.4577216) q[2];
sx q[2];
rz(-0.46087056) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8991124) q[1];
sx q[1];
rz(-0.63226262) q[1];
sx q[1];
rz(2.2805023) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71367587) q[3];
sx q[3];
rz(-1.6380966) q[3];
sx q[3];
rz(2.0590559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7632873) q[2];
sx q[2];
rz(-1.3102691) q[2];
sx q[2];
rz(2.7491167) q[2];
rz(1.1522419) q[3];
sx q[3];
rz(-2.4270054) q[3];
sx q[3];
rz(-0.31744441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1257989) q[0];
sx q[0];
rz(-1.5725461) q[0];
sx q[0];
rz(-2.3902067) q[0];
rz(-1.8136576) q[1];
sx q[1];
rz(-1.2633879) q[1];
sx q[1];
rz(0.60633916) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060121814) q[0];
sx q[0];
rz(-1.6166286) q[0];
sx q[0];
rz(1.5874552) q[0];
rz(-pi) q[1];
rz(2.2491127) q[2];
sx q[2];
rz(-1.2391029) q[2];
sx q[2];
rz(1.734037) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.025758) q[1];
sx q[1];
rz(-2.5587213) q[1];
sx q[1];
rz(1.9028266) q[1];
x q[2];
rz(-0.074999853) q[3];
sx q[3];
rz(-1.3538176) q[3];
sx q[3];
rz(2.0892339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5027344) q[2];
sx q[2];
rz(-2.0998462) q[2];
sx q[2];
rz(1.139337) q[2];
rz(-1.6566488) q[3];
sx q[3];
rz(-1.1805725) q[3];
sx q[3];
rz(-0.10425723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5320324) q[0];
sx q[0];
rz(-2.39344) q[0];
sx q[0];
rz(2.6334921) q[0];
rz(-1.5787026) q[1];
sx q[1];
rz(-2.0527614) q[1];
sx q[1];
rz(0.79024822) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9433141) q[0];
sx q[0];
rz(-2.4662848) q[0];
sx q[0];
rz(-0.4839464) q[0];
x q[1];
rz(3.087567) q[2];
sx q[2];
rz(-2.9276491) q[2];
sx q[2];
rz(0.33188785) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8040647) q[1];
sx q[1];
rz(-1.84066) q[1];
sx q[1];
rz(0.40562628) q[1];
rz(-pi) q[2];
rz(0.4831794) q[3];
sx q[3];
rz(-0.70671591) q[3];
sx q[3];
rz(-0.28373517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0885075) q[2];
sx q[2];
rz(-0.44184703) q[2];
sx q[2];
rz(-1.7283758) q[2];
rz(1.4767856) q[3];
sx q[3];
rz(-1.0352742) q[3];
sx q[3];
rz(-3.0800381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7168032) q[0];
sx q[0];
rz(-0.030310832) q[0];
sx q[0];
rz(-2.0943663) q[0];
rz(-0.60910243) q[1];
sx q[1];
rz(-1.7276238) q[1];
sx q[1];
rz(-1.3887127) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0213288) q[0];
sx q[0];
rz(-1.6006032) q[0];
sx q[0];
rz(2.1304312) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65865626) q[2];
sx q[2];
rz(-2.5052862) q[2];
sx q[2];
rz(-0.051740019) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.74275201) q[1];
sx q[1];
rz(-1.5821777) q[1];
sx q[1];
rz(-3.0497754) q[1];
x q[2];
rz(2.8453313) q[3];
sx q[3];
rz(-2.1564266) q[3];
sx q[3];
rz(-2.0010009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9528815) q[2];
sx q[2];
rz(-0.41023508) q[2];
sx q[2];
rz(0.88225538) q[2];
rz(1.4011718) q[3];
sx q[3];
rz(-1.975235) q[3];
sx q[3];
rz(-1.9410979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27154487) q[0];
sx q[0];
rz(-0.4168059) q[0];
sx q[0];
rz(1.7154988) q[0];
rz(3.0601314) q[1];
sx q[1];
rz(-1.9790244) q[1];
sx q[1];
rz(0.55823278) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0840069) q[0];
sx q[0];
rz(-1.1450197) q[0];
sx q[0];
rz(-2.7140679) q[0];
rz(1.4633281) q[2];
sx q[2];
rz(-1.5170013) q[2];
sx q[2];
rz(1.4449643) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0135632) q[1];
sx q[1];
rz(-0.82847825) q[1];
sx q[1];
rz(0.086704266) q[1];
x q[2];
rz(0.48645143) q[3];
sx q[3];
rz(-1.7351741) q[3];
sx q[3];
rz(0.96729507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2925064) q[2];
sx q[2];
rz(-1.2693274) q[2];
sx q[2];
rz(0.212184) q[2];
rz(-0.21197453) q[3];
sx q[3];
rz(-0.68325716) q[3];
sx q[3];
rz(1.9395444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50487173) q[0];
sx q[0];
rz(-0.8240521) q[0];
sx q[0];
rz(1.6037534) q[0];
rz(-2.3161855) q[1];
sx q[1];
rz(-0.67276612) q[1];
sx q[1];
rz(-0.5232946) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16738811) q[0];
sx q[0];
rz(-0.61199576) q[0];
sx q[0];
rz(0.1083072) q[0];
rz(-pi) q[1];
rz(1.5027572) q[2];
sx q[2];
rz(-0.57300742) q[2];
sx q[2];
rz(-2.9266561) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.20269468) q[1];
sx q[1];
rz(-1.885186) q[1];
sx q[1];
rz(2.7640192) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8966122) q[3];
sx q[3];
rz(-1.7953201) q[3];
sx q[3];
rz(-2.6905439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3045197) q[2];
sx q[2];
rz(-0.7154811) q[2];
sx q[2];
rz(-0.26930299) q[2];
rz(-2.6473911) q[3];
sx q[3];
rz(-0.84635693) q[3];
sx q[3];
rz(2.0555029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.3257278) q[0];
sx q[0];
rz(-1.5300735) q[0];
sx q[0];
rz(-1.6515401) q[0];
rz(1.6745463) q[1];
sx q[1];
rz(-0.29232262) q[1];
sx q[1];
rz(1.2437337) q[1];
rz(1.833563) q[2];
sx q[2];
rz(-1.5599712) q[2];
sx q[2];
rz(0.47906265) q[2];
rz(0.79694637) q[3];
sx q[3];
rz(-1.4016101) q[3];
sx q[3];
rz(-2.3102643) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];