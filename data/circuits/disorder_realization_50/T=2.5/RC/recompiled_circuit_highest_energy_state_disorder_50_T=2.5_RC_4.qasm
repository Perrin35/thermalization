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
rz(-1.4122352) q[0];
sx q[0];
rz(-0.5923624) q[0];
sx q[0];
rz(1.9368197) q[0];
rz(-2.0545948) q[1];
sx q[1];
rz(-0.49538651) q[1];
sx q[1];
rz(1.9538716) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43150768) q[0];
sx q[0];
rz(-2.1174701) q[0];
sx q[0];
rz(1.7492848) q[0];
x q[1];
rz(-2.4815791) q[2];
sx q[2];
rz(-0.68049351) q[2];
sx q[2];
rz(1.4726435) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4036619) q[1];
sx q[1];
rz(-1.523535) q[1];
sx q[1];
rz(-1.74356) q[1];
x q[2];
rz(1.2212444) q[3];
sx q[3];
rz(-0.38844019) q[3];
sx q[3];
rz(0.74439186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5997233) q[2];
sx q[2];
rz(-2.5399127) q[2];
sx q[2];
rz(2.0987233) q[2];
rz(-3.0964105) q[3];
sx q[3];
rz(-2.9729645) q[3];
sx q[3];
rz(-1.5359623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0586108) q[0];
sx q[0];
rz(-3.0284212) q[0];
sx q[0];
rz(2.163072) q[0];
rz(1.9784031) q[1];
sx q[1];
rz(-0.99449831) q[1];
sx q[1];
rz(1.3226343) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.037905692) q[0];
sx q[0];
rz(-0.9390489) q[0];
sx q[0];
rz(-0.61236713) q[0];
rz(-pi) q[1];
rz(2.1345226) q[2];
sx q[2];
rz(-0.33180922) q[2];
sx q[2];
rz(-0.72124583) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1910227) q[1];
sx q[1];
rz(-0.92838597) q[1];
sx q[1];
rz(2.8710142) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7197116) q[3];
sx q[3];
rz(-2.2595539) q[3];
sx q[3];
rz(1.0721579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.80403745) q[2];
sx q[2];
rz(-0.20122169) q[2];
sx q[2];
rz(0.03841722) q[2];
rz(-1.6328968) q[3];
sx q[3];
rz(-1.8149866) q[3];
sx q[3];
rz(1.4940777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22817336) q[0];
sx q[0];
rz(-2.4689624) q[0];
sx q[0];
rz(0.2359373) q[0];
rz(-1.1783696) q[1];
sx q[1];
rz(-1.6940247) q[1];
sx q[1];
rz(-0.49740121) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.936393) q[0];
sx q[0];
rz(-2.8872364) q[0];
sx q[0];
rz(0.78820552) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3942986) q[2];
sx q[2];
rz(-1.1836393) q[2];
sx q[2];
rz(2.5839069) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6420501) q[1];
sx q[1];
rz(-1.9643557) q[1];
sx q[1];
rz(-1.5469013) q[1];
rz(-pi) q[2];
x q[2];
rz(2.892832) q[3];
sx q[3];
rz(-1.3205055) q[3];
sx q[3];
rz(0.32969013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.95911372) q[2];
sx q[2];
rz(-1.6751869) q[2];
sx q[2];
rz(1.2314931) q[2];
rz(1.4500827) q[3];
sx q[3];
rz(-1.8385889) q[3];
sx q[3];
rz(0.85795295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0602144) q[0];
sx q[0];
rz(-0.84351081) q[0];
sx q[0];
rz(2.9486935) q[0];
rz(1.3554205) q[1];
sx q[1];
rz(-1.5813634) q[1];
sx q[1];
rz(-0.17098175) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9381704) q[0];
sx q[0];
rz(-1.6294384) q[0];
sx q[0];
rz(-1.0542271) q[0];
x q[1];
rz(1.0320206) q[2];
sx q[2];
rz(-2.2316885) q[2];
sx q[2];
rz(0.62884841) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0031922) q[1];
sx q[1];
rz(-2.698771) q[1];
sx q[1];
rz(-1.3243616) q[1];
rz(-pi) q[2];
rz(0.2410197) q[3];
sx q[3];
rz(-2.2200826) q[3];
sx q[3];
rz(-2.5229682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1811447) q[2];
sx q[2];
rz(-1.9713216) q[2];
sx q[2];
rz(2.81874) q[2];
rz(1.3747831) q[3];
sx q[3];
rz(-2.1026473) q[3];
sx q[3];
rz(2.3006181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5508995) q[0];
sx q[0];
rz(-2.1192079) q[0];
sx q[0];
rz(-0.92556959) q[0];
rz(-3.0135221) q[1];
sx q[1];
rz(-1.4984727) q[1];
sx q[1];
rz(-3.0421323) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8443237) q[0];
sx q[0];
rz(-0.8542866) q[0];
sx q[0];
rz(-3.141527) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2242618) q[2];
sx q[2];
rz(-1.230403) q[2];
sx q[2];
rz(0.16763359) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4872553) q[1];
sx q[1];
rz(-1.0434101) q[1];
sx q[1];
rz(2.9156963) q[1];
x q[2];
rz(-3.0614047) q[3];
sx q[3];
rz(-1.8155451) q[3];
sx q[3];
rz(-0.65746869) q[3];
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
rz(-1.1348178) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9099092) q[0];
sx q[0];
rz(-2.9250356) q[0];
sx q[0];
rz(-0.59189558) q[0];
rz(-0.20467155) q[1];
sx q[1];
rz(-0.99212956) q[1];
sx q[1];
rz(3.0904904) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0436193) q[0];
sx q[0];
rz(-1.4187001) q[0];
sx q[0];
rz(-0.19466227) q[0];
rz(-pi) q[1];
rz(-2.2289071) q[2];
sx q[2];
rz(-1.776623) q[2];
sx q[2];
rz(1.4806946) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.49858677) q[1];
sx q[1];
rz(-2.6664147) q[1];
sx q[1];
rz(-2.2982486) q[1];
rz(-pi) q[2];
x q[2];
rz(0.74441142) q[3];
sx q[3];
rz(-2.4648414) q[3];
sx q[3];
rz(1.6729421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0245612) q[2];
sx q[2];
rz(-2.0270963) q[2];
sx q[2];
rz(2.047211) q[2];
rz(-0.38703212) q[3];
sx q[3];
rz(-2.6670167) q[3];
sx q[3];
rz(-0.3064557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36560202) q[0];
sx q[0];
rz(-2.2633573) q[0];
sx q[0];
rz(0.24755092) q[0];
rz(1.4546825) q[1];
sx q[1];
rz(-1.8984112) q[1];
sx q[1];
rz(-3.1051342) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.582583) q[0];
sx q[0];
rz(-1.2866308) q[0];
sx q[0];
rz(-0.76822944) q[0];
rz(-pi) q[1];
rz(1.7629303) q[2];
sx q[2];
rz(-2.2647144) q[2];
sx q[2];
rz(-1.0200227) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.16462886) q[1];
sx q[1];
rz(-2.0348624) q[1];
sx q[1];
rz(1.0010757) q[1];
x q[2];
rz(3.0667449) q[3];
sx q[3];
rz(-1.6597431) q[3];
sx q[3];
rz(-1.4527827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0869007) q[2];
sx q[2];
rz(-0.96618235) q[2];
sx q[2];
rz(-1.7894233) q[2];
rz(-1.8933206) q[3];
sx q[3];
rz(-1.5779481) q[3];
sx q[3];
rz(-1.9498391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9098814) q[0];
sx q[0];
rz(-0.24869643) q[0];
sx q[0];
rz(1.4924208) q[0];
rz(2.9630648) q[1];
sx q[1];
rz(-1.2622204) q[1];
sx q[1];
rz(2.5915204) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72075413) q[0];
sx q[0];
rz(-1.5394569) q[0];
sx q[0];
rz(3.1216338) q[0];
rz(-1.3230611) q[2];
sx q[2];
rz(-1.6322517) q[2];
sx q[2];
rz(2.2814192) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7272721) q[1];
sx q[1];
rz(-1.7309233) q[1];
sx q[1];
rz(-0.54267197) q[1];
rz(-pi) q[2];
rz(-2.7356803) q[3];
sx q[3];
rz(-1.1391057) q[3];
sx q[3];
rz(-1.622705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7940346) q[2];
sx q[2];
rz(-1.4057691) q[2];
sx q[2];
rz(-2.3883635) q[2];
rz(-2.8473162) q[3];
sx q[3];
rz(-0.080065057) q[3];
sx q[3];
rz(0.068597138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8681965) q[0];
sx q[0];
rz(-0.28565872) q[0];
sx q[0];
rz(-1.4519325) q[0];
rz(-0.1952576) q[1];
sx q[1];
rz(-2.0611019) q[1];
sx q[1];
rz(1.4917699) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0367835) q[0];
sx q[0];
rz(-1.097479) q[0];
sx q[0];
rz(-1.2261054) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3765747) q[2];
sx q[2];
rz(-1.6510626) q[2];
sx q[2];
rz(1.8945872) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0635707) q[1];
sx q[1];
rz(-1.6352904) q[1];
sx q[1];
rz(1.4149354) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5641065) q[3];
sx q[3];
rz(-1.5433528) q[3];
sx q[3];
rz(0.60718482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1008272) q[2];
sx q[2];
rz(-2.6746174) q[2];
sx q[2];
rz(-0.85774285) q[2];
rz(1.6929251) q[3];
sx q[3];
rz(-1.6489776) q[3];
sx q[3];
rz(2.9878476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.0070852) q[0];
sx q[0];
rz(-0.15468287) q[0];
sx q[0];
rz(2.3585228) q[0];
rz(2.018441) q[1];
sx q[1];
rz(-1.4479366) q[1];
sx q[1];
rz(3.0812841) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9674112) q[0];
sx q[0];
rz(-2.7740712) q[0];
sx q[0];
rz(1.307748) q[0];
rz(-1.7293634) q[2];
sx q[2];
rz(-2.5494908) q[2];
sx q[2];
rz(-0.71449661) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.26877182) q[1];
sx q[1];
rz(-0.67292456) q[1];
sx q[1];
rz(-3.0329143) q[1];
rz(-pi) q[2];
rz(2.7785886) q[3];
sx q[3];
rz(-1.6771183) q[3];
sx q[3];
rz(-2.8385988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15410885) q[2];
sx q[2];
rz(-2.2769112) q[2];
sx q[2];
rz(-1.8314499) q[2];
rz(1.3335258) q[3];
sx q[3];
rz(-1.798505) q[3];
sx q[3];
rz(-1.8346627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46351984) q[0];
sx q[0];
rz(-1.1165883) q[0];
sx q[0];
rz(-0.22620871) q[0];
rz(-0.75640596) q[1];
sx q[1];
rz(-0.49497985) q[1];
sx q[1];
rz(2.3351647) q[1];
rz(0.020515223) q[2];
sx q[2];
rz(-2.7111369) q[2];
sx q[2];
rz(-2.6028056) q[2];
rz(-1.9202833) q[3];
sx q[3];
rz(-1.7428453) q[3];
sx q[3];
rz(1.3475628) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
