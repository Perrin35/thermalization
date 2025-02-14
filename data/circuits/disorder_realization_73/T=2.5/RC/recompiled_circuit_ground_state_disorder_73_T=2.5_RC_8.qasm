OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.29326987) q[0];
sx q[0];
rz(-2.8889416) q[0];
sx q[0];
rz(1.2055612) q[0];
rz(2.0134917) q[1];
sx q[1];
rz(-1.3265346) q[1];
sx q[1];
rz(0.74572745) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0416314) q[0];
sx q[0];
rz(-1.7929165) q[0];
sx q[0];
rz(-3.0801386) q[0];
rz(-pi) q[1];
rz(-3.1118561) q[2];
sx q[2];
rz(-0.22448891) q[2];
sx q[2];
rz(-2.9475074) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9348976) q[1];
sx q[1];
rz(-1.3612011) q[1];
sx q[1];
rz(1.2278201) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9609591) q[3];
sx q[3];
rz(-2.6350065) q[3];
sx q[3];
rz(2.6449639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.91743177) q[2];
sx q[2];
rz(-1.5956343) q[2];
sx q[2];
rz(1.6708299) q[2];
rz(-1.6338232) q[3];
sx q[3];
rz(-3.1247415) q[3];
sx q[3];
rz(0.9233709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5725752) q[0];
sx q[0];
rz(-1.1976396) q[0];
sx q[0];
rz(-1.5684599) q[0];
rz(-2.9720427) q[1];
sx q[1];
rz(-3.0270271) q[1];
sx q[1];
rz(-0.13470185) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3189663) q[0];
sx q[0];
rz(-1.14577) q[0];
sx q[0];
rz(-1.1898196) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6394272) q[2];
sx q[2];
rz(-1.5585715) q[2];
sx q[2];
rz(3.1032012) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0162279) q[1];
sx q[1];
rz(-1.1332336) q[1];
sx q[1];
rz(-1.2368278) q[1];
x q[2];
rz(-1.991113) q[3];
sx q[3];
rz(-0.1945217) q[3];
sx q[3];
rz(1.9357301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0160825) q[2];
sx q[2];
rz(-1.4849911) q[2];
sx q[2];
rz(-2.9888195) q[2];
rz(-1.3656535) q[3];
sx q[3];
rz(-3.1047265) q[3];
sx q[3];
rz(-2.9735145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.633054) q[0];
sx q[0];
rz(-0.80798739) q[0];
sx q[0];
rz(0.48164865) q[0];
rz(-2.956849) q[1];
sx q[1];
rz(-1.3667204) q[1];
sx q[1];
rz(0.99536037) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5115968) q[0];
sx q[0];
rz(-2.7212486) q[0];
sx q[0];
rz(-2.1612801) q[0];
rz(-0.89389385) q[2];
sx q[2];
rz(-3.0640996) q[2];
sx q[2];
rz(-0.95439664) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3832356) q[1];
sx q[1];
rz(-2.2192973) q[1];
sx q[1];
rz(1.1926257) q[1];
x q[2];
rz(1.6158726) q[3];
sx q[3];
rz(-0.59642506) q[3];
sx q[3];
rz(0.10208043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92622009) q[2];
sx q[2];
rz(-0.054457713) q[2];
sx q[2];
rz(-0.064662956) q[2];
rz(2.0926545) q[3];
sx q[3];
rz(-0.026853042) q[3];
sx q[3];
rz(-1.8612727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8365086) q[0];
sx q[0];
rz(-3.0236112) q[0];
sx q[0];
rz(-0.82103658) q[0];
rz(-3.0631284) q[1];
sx q[1];
rz(-1.6946225) q[1];
sx q[1];
rz(-2.1441114) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8686912) q[0];
sx q[0];
rz(-1.8857546) q[0];
sx q[0];
rz(0.64906831) q[0];
x q[1];
rz(1.5115963) q[2];
sx q[2];
rz(-1.6069876) q[2];
sx q[2];
rz(0.22360392) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9033244) q[1];
sx q[1];
rz(-2.0966171) q[1];
sx q[1];
rz(2.8366798) q[1];
rz(-3.0313052) q[3];
sx q[3];
rz(-0.48088851) q[3];
sx q[3];
rz(2.4824924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0923826) q[2];
sx q[2];
rz(-0.02928484) q[2];
sx q[2];
rz(-1.9878261) q[2];
rz(2.9554101) q[3];
sx q[3];
rz(-3.0598873) q[3];
sx q[3];
rz(2.8103099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.004772923) q[0];
sx q[0];
rz(-0.81757075) q[0];
sx q[0];
rz(-1.1432884) q[0];
rz(2.000287) q[1];
sx q[1];
rz(-2.3316796) q[1];
sx q[1];
rz(0.51796651) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2835613) q[0];
sx q[0];
rz(-2.1037344) q[0];
sx q[0];
rz(-2.3951247) q[0];
x q[1];
rz(-1.5094366) q[2];
sx q[2];
rz(-3.1265335) q[2];
sx q[2];
rz(1.7565774) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.28455602) q[1];
sx q[1];
rz(-1.9003881) q[1];
sx q[1];
rz(-0.28917851) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7180743) q[3];
sx q[3];
rz(-1.8657576) q[3];
sx q[3];
rz(-1.1671821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3190069) q[2];
sx q[2];
rz(-0.95390445) q[2];
sx q[2];
rz(0.32773584) q[2];
rz(-1.1945126) q[3];
sx q[3];
rz(-0.12735282) q[3];
sx q[3];
rz(-2.2849042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1389393) q[0];
sx q[0];
rz(-2.8728573) q[0];
sx q[0];
rz(0.56185454) q[0];
rz(1.6506763) q[1];
sx q[1];
rz(-1.6249388) q[1];
sx q[1];
rz(-3.0432826) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5033103) q[0];
sx q[0];
rz(-1.1489963) q[0];
sx q[0];
rz(0.038859239) q[0];
rz(-pi) q[1];
rz(-0.00064050015) q[2];
sx q[2];
rz(-1.5709236) q[2];
sx q[2];
rz(-1.887111) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1001818) q[1];
sx q[1];
rz(-1.6173925) q[1];
sx q[1];
rz(-0.46943922) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6410511) q[3];
sx q[3];
rz(-1.9841474) q[3];
sx q[3];
rz(-2.7374637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.085122846) q[2];
sx q[2];
rz(-2.9462908) q[2];
sx q[2];
rz(-1.0657715) q[2];
rz(-0.39984518) q[3];
sx q[3];
rz(-2.6078434) q[3];
sx q[3];
rz(-1.2679509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(0.14168508) q[0];
sx q[0];
rz(-3.0596924) q[0];
sx q[0];
rz(1.4578693) q[0];
rz(1.1204002) q[1];
sx q[1];
rz(-3.0034062) q[1];
sx q[1];
rz(-2.8021326) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0908112) q[0];
sx q[0];
rz(-3.0635186) q[0];
sx q[0];
rz(-0.17103057) q[0];
rz(-pi) q[1];
rz(-3.1192084) q[2];
sx q[2];
rz(-0.55479344) q[2];
sx q[2];
rz(-3.1172397) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5684699) q[1];
sx q[1];
rz(-3.0559799) q[1];
sx q[1];
rz(-1.5898311) q[1];
x q[2];
rz(2.8651627) q[3];
sx q[3];
rz(-2.7904841) q[3];
sx q[3];
rz(2.1123667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7772943) q[2];
sx q[2];
rz(-3.1396781) q[2];
sx q[2];
rz(-0.36336362) q[2];
rz(2.0511138) q[3];
sx q[3];
rz(-2.5616779) q[3];
sx q[3];
rz(-1.1183848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5217487) q[0];
sx q[0];
rz(-0.32829568) q[0];
sx q[0];
rz(2.2186665) q[0];
rz(1.482831) q[1];
sx q[1];
rz(-2.5117579) q[1];
sx q[1];
rz(-3.1373851) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.370458) q[0];
sx q[0];
rz(-2.1677289) q[0];
sx q[0];
rz(-2.202522) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5571874) q[2];
sx q[2];
rz(-0.41092349) q[2];
sx q[2];
rz(3.1142554) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1394704) q[1];
sx q[1];
rz(-1.0957901) q[1];
sx q[1];
rz(-0.072474555) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.065807314) q[3];
sx q[3];
rz(-2.5281457) q[3];
sx q[3];
rz(-2.891401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5620455) q[2];
sx q[2];
rz(-1.5374708) q[2];
sx q[2];
rz(1.9407678) q[2];
rz(1.3935401) q[3];
sx q[3];
rz(-3.1378919) q[3];
sx q[3];
rz(0.70011955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8414108) q[0];
sx q[0];
rz(-0.44898471) q[0];
sx q[0];
rz(-1.9218943) q[0];
rz(-1.8118743) q[1];
sx q[1];
rz(-2.0025573) q[1];
sx q[1];
rz(-2.9776998) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3541479) q[0];
sx q[0];
rz(-1.3992926) q[0];
sx q[0];
rz(-1.9308912) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40028769) q[2];
sx q[2];
rz(-2.5185555) q[2];
sx q[2];
rz(-1.7965339) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.51162219) q[1];
sx q[1];
rz(-3.0093806) q[1];
sx q[1];
rz(-1.1019143) q[1];
rz(-pi) q[2];
rz(-2.4285942) q[3];
sx q[3];
rz(-2.2621691) q[3];
sx q[3];
rz(-0.52934968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4280052) q[2];
sx q[2];
rz(-0.089944936) q[2];
sx q[2];
rz(0.62291992) q[2];
rz(-0.04341393) q[3];
sx q[3];
rz(-2.2446938) q[3];
sx q[3];
rz(-2.3626732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49122214) q[0];
sx q[0];
rz(-3.1351884) q[0];
sx q[0];
rz(-2.6553335) q[0];
rz(-2.4468415) q[1];
sx q[1];
rz(-0.34725747) q[1];
sx q[1];
rz(-0.4756701) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2025288) q[0];
sx q[0];
rz(-0.049204218) q[0];
sx q[0];
rz(2.3243107) q[0];
rz(-pi) q[1];
rz(2.2769502) q[2];
sx q[2];
rz(-1.5967895) q[2];
sx q[2];
rz(-0.038250462) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5359395) q[1];
sx q[1];
rz(-1.5772784) q[1];
sx q[1];
rz(-2.5340213) q[1];
rz(-pi) q[2];
rz(0.31020152) q[3];
sx q[3];
rz(-1.3942914) q[3];
sx q[3];
rz(-1.4909397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4645369) q[2];
sx q[2];
rz(-3.0812283) q[2];
sx q[2];
rz(1.2976868) q[2];
rz(-1.8929947) q[3];
sx q[3];
rz(-0.55515754) q[3];
sx q[3];
rz(-2.7738074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1260592) q[0];
sx q[0];
rz(-1.5801237) q[0];
sx q[0];
rz(1.7042241) q[0];
rz(-2.2592648) q[1];
sx q[1];
rz(-3.075141) q[1];
sx q[1];
rz(-2.3022423) q[1];
rz(-2.7693314) q[2];
sx q[2];
rz(-3.0426171) q[2];
sx q[2];
rz(-0.097886861) q[2];
rz(-1.5658548) q[3];
sx q[3];
rz(-1.3316122) q[3];
sx q[3];
rz(0.022461654) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
