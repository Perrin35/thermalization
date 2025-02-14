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
rz(-2.8430804) q[0];
sx q[0];
rz(-1.2895583) q[0];
sx q[0];
rz(-0.02331743) q[0];
rz(-2.159637) q[1];
sx q[1];
rz(-2.4040931) q[1];
sx q[1];
rz(-1.4758543) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1630572) q[0];
sx q[0];
rz(-0.082162372) q[0];
sx q[0];
rz(2.1294247) q[0];
x q[1];
rz(1.2405922) q[2];
sx q[2];
rz(-1.6088952) q[2];
sx q[2];
rz(2.7862501) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.45881328) q[1];
sx q[1];
rz(-0.87535697) q[1];
sx q[1];
rz(1.0473482) q[1];
rz(-pi) q[2];
rz(-0.89453434) q[3];
sx q[3];
rz(-1.1522376) q[3];
sx q[3];
rz(0.422131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.6797356) q[2];
sx q[2];
rz(-1.6697786) q[2];
sx q[2];
rz(0.30065817) q[2];
rz(2.4833208) q[3];
sx q[3];
rz(-2.7418147) q[3];
sx q[3];
rz(-2.5532653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7268426) q[0];
sx q[0];
rz(-2.2756133) q[0];
sx q[0];
rz(-0.17587371) q[0];
rz(2.580592) q[1];
sx q[1];
rz(-2.439552) q[1];
sx q[1];
rz(-2.9452513) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8844159) q[0];
sx q[0];
rz(-1.1286583) q[0];
sx q[0];
rz(-0.37143645) q[0];
x q[1];
rz(-3.0220993) q[2];
sx q[2];
rz(-1.4548848) q[2];
sx q[2];
rz(1.3136065) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.93028101) q[1];
sx q[1];
rz(-2.086211) q[1];
sx q[1];
rz(-0.0023018259) q[1];
x q[2];
rz(1.4370055) q[3];
sx q[3];
rz(-1.2302629) q[3];
sx q[3];
rz(1.2029369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0520797) q[2];
sx q[2];
rz(-0.61605805) q[2];
sx q[2];
rz(1.6717795) q[2];
rz(0.52516627) q[3];
sx q[3];
rz(-1.7142121) q[3];
sx q[3];
rz(-2.4770881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1216275) q[0];
sx q[0];
rz(-0.88931924) q[0];
sx q[0];
rz(2.5110733) q[0];
rz(1.9969253) q[1];
sx q[1];
rz(-1.8960709) q[1];
sx q[1];
rz(-3.0847881) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58724126) q[0];
sx q[0];
rz(-1.2040753) q[0];
sx q[0];
rz(1.1318867) q[0];
x q[1];
rz(0.47614758) q[2];
sx q[2];
rz(-2.2000929) q[2];
sx q[2];
rz(-0.72222939) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.60384974) q[1];
sx q[1];
rz(-1.4293539) q[1];
sx q[1];
rz(-0.46700041) q[1];
rz(-pi) q[2];
rz(0.2127368) q[3];
sx q[3];
rz(-1.9895344) q[3];
sx q[3];
rz(-1.9902844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8593665) q[2];
sx q[2];
rz(-1.4713918) q[2];
sx q[2];
rz(-2.7024506) q[2];
rz(1.7450843) q[3];
sx q[3];
rz(-1.8789995) q[3];
sx q[3];
rz(0.5784353) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0466995) q[0];
sx q[0];
rz(-2.4201604) q[0];
sx q[0];
rz(2.1266345) q[0];
rz(-0.062648423) q[1];
sx q[1];
rz(-0.95476127) q[1];
sx q[1];
rz(0.28422022) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64062041) q[0];
sx q[0];
rz(-0.56864244) q[0];
sx q[0];
rz(0.68822648) q[0];
rz(-pi) q[1];
rz(-2.6437628) q[2];
sx q[2];
rz(-1.6938899) q[2];
sx q[2];
rz(0.13409889) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9584459) q[1];
sx q[1];
rz(-1.5827521) q[1];
sx q[1];
rz(2.7832116) q[1];
rz(-pi) q[2];
rz(-1.5169859) q[3];
sx q[3];
rz(-2.887886) q[3];
sx q[3];
rz(1.0347317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8404428) q[2];
sx q[2];
rz(-1.0089077) q[2];
sx q[2];
rz(-0.76048771) q[2];
rz(2.8216951) q[3];
sx q[3];
rz(-1.1287929) q[3];
sx q[3];
rz(-2.1130051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8496721) q[0];
sx q[0];
rz(-0.55067647) q[0];
sx q[0];
rz(0.40529761) q[0];
rz(0.48794508) q[1];
sx q[1];
rz(-2.0351724) q[1];
sx q[1];
rz(-0.36283666) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3170119) q[0];
sx q[0];
rz(-1.2634227) q[0];
sx q[0];
rz(2.6885343) q[0];
rz(-pi) q[1];
rz(-1.9486729) q[2];
sx q[2];
rz(-2.123017) q[2];
sx q[2];
rz(-0.044980031) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.47687864) q[1];
sx q[1];
rz(-1.6194317) q[1];
sx q[1];
rz(-0.22108476) q[1];
rz(0.36146253) q[3];
sx q[3];
rz(-0.72849724) q[3];
sx q[3];
rz(-1.1379153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1399347) q[2];
sx q[2];
rz(-1.9826865) q[2];
sx q[2];
rz(-0.72251594) q[2];
rz(2.6288988) q[3];
sx q[3];
rz(-0.86944681) q[3];
sx q[3];
rz(-1.2041913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8743415) q[0];
sx q[0];
rz(-2.5047472) q[0];
sx q[0];
rz(0.27477086) q[0];
rz(0.35708669) q[1];
sx q[1];
rz(-1.6386702) q[1];
sx q[1];
rz(1.1579317) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5362893) q[0];
sx q[0];
rz(-2.2380296) q[0];
sx q[0];
rz(0.31179223) q[0];
rz(-pi) q[1];
rz(-2.2311796) q[2];
sx q[2];
rz(-2.0636301) q[2];
sx q[2];
rz(3.1255258) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.14640644) q[1];
sx q[1];
rz(-2.1651092) q[1];
sx q[1];
rz(1.8176778) q[1];
x q[2];
rz(2.3401392) q[3];
sx q[3];
rz(-1.0545925) q[3];
sx q[3];
rz(-3.0954697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5623515) q[2];
sx q[2];
rz(-1.7513821) q[2];
sx q[2];
rz(-0.58970279) q[2];
rz(1.8288745) q[3];
sx q[3];
rz(-1.6629985) q[3];
sx q[3];
rz(0.5635128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24484672) q[0];
sx q[0];
rz(-0.99269301) q[0];
sx q[0];
rz(-2.8712811) q[0];
rz(-0.42701328) q[1];
sx q[1];
rz(-1.7662363) q[1];
sx q[1];
rz(0.40831533) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9664551) q[0];
sx q[0];
rz(-1.2027367) q[0];
sx q[0];
rz(1.6987263) q[0];
rz(3.0403942) q[2];
sx q[2];
rz(-2.2658341) q[2];
sx q[2];
rz(1.1388495) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4998251) q[1];
sx q[1];
rz(-2.1534792) q[1];
sx q[1];
rz(0.50104435) q[1];
rz(-pi) q[2];
rz(1.6736843) q[3];
sx q[3];
rz(-0.85082084) q[3];
sx q[3];
rz(0.83482367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2030877) q[2];
sx q[2];
rz(-2.0577343) q[2];
sx q[2];
rz(-2.2173524) q[2];
rz(0.83485323) q[3];
sx q[3];
rz(-2.3122841) q[3];
sx q[3];
rz(-0.99669641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95558178) q[0];
sx q[0];
rz(-2.7614433) q[0];
sx q[0];
rz(2.7741449) q[0];
rz(2.5732749) q[1];
sx q[1];
rz(-1.3325997) q[1];
sx q[1];
rz(-2.7265991) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5552705) q[0];
sx q[0];
rz(-1.30997) q[0];
sx q[0];
rz(-1.1010948) q[0];
rz(2.1176362) q[2];
sx q[2];
rz(-2.0342361) q[2];
sx q[2];
rz(-1.5228423) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3967949) q[1];
sx q[1];
rz(-0.84229453) q[1];
sx q[1];
rz(-0.70835967) q[1];
x q[2];
rz(-1.2250118) q[3];
sx q[3];
rz(-1.5659075) q[3];
sx q[3];
rz(1.1556311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.35855287) q[2];
sx q[2];
rz(-0.8873322) q[2];
sx q[2];
rz(-2.3084194) q[2];
rz(-2.7821275) q[3];
sx q[3];
rz(-1.0901674) q[3];
sx q[3];
rz(1.5045213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-1.9678765) q[0];
sx q[0];
rz(-0.40206566) q[0];
sx q[0];
rz(0.014884431) q[0];
rz(-2.0170276) q[1];
sx q[1];
rz(-2.896307) q[1];
sx q[1];
rz(3.0344149) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2994605) q[0];
sx q[0];
rz(-0.74343527) q[0];
sx q[0];
rz(1.7920262) q[0];
x q[1];
rz(-1.0753267) q[2];
sx q[2];
rz(-0.81186724) q[2];
sx q[2];
rz(2.5798485) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5276423) q[1];
sx q[1];
rz(-2.6384934) q[1];
sx q[1];
rz(-0.31219074) q[1];
rz(1.7483057) q[3];
sx q[3];
rz(-1.313398) q[3];
sx q[3];
rz(0.37198113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5643481) q[2];
sx q[2];
rz(-1.3078835) q[2];
sx q[2];
rz(0.85961071) q[2];
rz(2.882242) q[3];
sx q[3];
rz(-1.3602463) q[3];
sx q[3];
rz(-0.41659659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49941007) q[0];
sx q[0];
rz(-3.0152617) q[0];
sx q[0];
rz(-1.9835749) q[0];
rz(-2.6462789) q[1];
sx q[1];
rz(-1.1664349) q[1];
sx q[1];
rz(-2.8445629) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3756145) q[0];
sx q[0];
rz(-2.8225464) q[0];
sx q[0];
rz(1.5578397) q[0];
rz(-pi) q[1];
rz(1.530324) q[2];
sx q[2];
rz(-1.8708611) q[2];
sx q[2];
rz(1.7057246) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90290711) q[1];
sx q[1];
rz(-0.33329646) q[1];
sx q[1];
rz(2.6413647) q[1];
x q[2];
rz(-0.90822753) q[3];
sx q[3];
rz(-1.5836827) q[3];
sx q[3];
rz(2.9290309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.26934066) q[2];
sx q[2];
rz(-0.38186914) q[2];
sx q[2];
rz(-0.86088172) q[2];
rz(-2.6779029) q[3];
sx q[3];
rz(-2.2380565) q[3];
sx q[3];
rz(1.1369107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9853482) q[0];
sx q[0];
rz(-1.5679659) q[0];
sx q[0];
rz(1.1563942) q[0];
rz(-1.3337749) q[1];
sx q[1];
rz(-1.540624) q[1];
sx q[1];
rz(-2.435138) q[1];
rz(-2.479066) q[2];
sx q[2];
rz(-2.0338197) q[2];
sx q[2];
rz(2.8070569) q[2];
rz(-0.49336962) q[3];
sx q[3];
rz(-1.9153638) q[3];
sx q[3];
rz(2.3079688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
