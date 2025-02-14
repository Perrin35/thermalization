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
rz(-1.2920657) q[0];
sx q[0];
rz(-2.3176365) q[0];
sx q[0];
rz(-2.2141461) q[0];
rz(-1.0499586) q[1];
sx q[1];
rz(-2.1702622) q[1];
sx q[1];
rz(0.8144905) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8626635) q[0];
sx q[0];
rz(-1.188108) q[0];
sx q[0];
rz(2.0576059) q[0];
rz(-2.522722) q[2];
sx q[2];
rz(-2.2437988) q[2];
sx q[2];
rz(-0.25970632) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8100639) q[1];
sx q[1];
rz(-2.0671856) q[1];
sx q[1];
rz(1.0277102) q[1];
rz(-pi) q[2];
rz(-0.85759576) q[3];
sx q[3];
rz(-2.2603545) q[3];
sx q[3];
rz(-1.5201207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3689975) q[2];
sx q[2];
rz(-1.3694222) q[2];
sx q[2];
rz(-0.43130809) q[2];
rz(1.6469693) q[3];
sx q[3];
rz(-1.5957811) q[3];
sx q[3];
rz(-2.3561884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.5539219) q[0];
sx q[0];
rz(-2.9449154) q[0];
sx q[0];
rz(1.824463) q[0];
rz(1.0493086) q[1];
sx q[1];
rz(-0.53593719) q[1];
sx q[1];
rz(2.5228693) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4878889) q[0];
sx q[0];
rz(-1.426322) q[0];
sx q[0];
rz(2.5313951) q[0];
rz(-0.40561442) q[2];
sx q[2];
rz(-2.2219293) q[2];
sx q[2];
rz(1.3319911) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4977469) q[1];
sx q[1];
rz(-0.37393716) q[1];
sx q[1];
rz(-1.9356329) q[1];
rz(-pi) q[2];
rz(-1.9369164) q[3];
sx q[3];
rz(-1.06166) q[3];
sx q[3];
rz(2.3222271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.416136) q[2];
sx q[2];
rz(-1.8014896) q[2];
sx q[2];
rz(-0.49187342) q[2];
rz(-1.5623931) q[3];
sx q[3];
rz(-1.4926566) q[3];
sx q[3];
rz(-2.5621342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014515011) q[0];
sx q[0];
rz(-2.0836232) q[0];
sx q[0];
rz(-0.28011093) q[0];
rz(0.07490553) q[1];
sx q[1];
rz(-0.43498755) q[1];
sx q[1];
rz(1.3704971) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89557392) q[0];
sx q[0];
rz(-1.7798806) q[0];
sx q[0];
rz(-0.81935291) q[0];
x q[1];
rz(0.63997322) q[2];
sx q[2];
rz(-1.0037862) q[2];
sx q[2];
rz(-2.8550573) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4968003) q[1];
sx q[1];
rz(-2.519517) q[1];
sx q[1];
rz(0.46228564) q[1];
rz(-pi) q[2];
rz(-3.0093569) q[3];
sx q[3];
rz(-1.5894798) q[3];
sx q[3];
rz(-1.736369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.28477731) q[2];
sx q[2];
rz(-1.4692401) q[2];
sx q[2];
rz(-2.3397297) q[2];
rz(-2.2297468) q[3];
sx q[3];
rz(-2.3175479) q[3];
sx q[3];
rz(2.1626933) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80482471) q[0];
sx q[0];
rz(-2.473859) q[0];
sx q[0];
rz(2.5208852) q[0];
rz(-0.2785109) q[1];
sx q[1];
rz(-1.6997489) q[1];
sx q[1];
rz(-1.0034358) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1708978) q[0];
sx q[0];
rz(-1.8866294) q[0];
sx q[0];
rz(1.210808) q[0];
rz(-pi) q[1];
x q[1];
rz(1.150028) q[2];
sx q[2];
rz(-1.506797) q[2];
sx q[2];
rz(-0.45129946) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4651854) q[1];
sx q[1];
rz(-2.8670951) q[1];
sx q[1];
rz(2.9390915) q[1];
rz(0.33903709) q[3];
sx q[3];
rz(-0.62612306) q[3];
sx q[3];
rz(2.1563185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3044546) q[2];
sx q[2];
rz(-1.1123603) q[2];
sx q[2];
rz(-3.0755074) q[2];
rz(1.1262013) q[3];
sx q[3];
rz(-2.3812713) q[3];
sx q[3];
rz(2.5662305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5657625) q[0];
sx q[0];
rz(-0.081871651) q[0];
sx q[0];
rz(0.78731147) q[0];
rz(-0.85834223) q[1];
sx q[1];
rz(-2.1635677) q[1];
sx q[1];
rz(-1.0816921) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72321945) q[0];
sx q[0];
rz(-2.0670509) q[0];
sx q[0];
rz(0.70696522) q[0];
rz(-pi) q[1];
rz(-2.0400286) q[2];
sx q[2];
rz(-1.9479556) q[2];
sx q[2];
rz(-2.9531329) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8906403) q[1];
sx q[1];
rz(-1.8554652) q[1];
sx q[1];
rz(-1.5773058) q[1];
rz(-0.55629913) q[3];
sx q[3];
rz(-0.93141684) q[3];
sx q[3];
rz(-2.104708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9786238) q[2];
sx q[2];
rz(-1.1498412) q[2];
sx q[2];
rz(-0.84257379) q[2];
rz(0.27635559) q[3];
sx q[3];
rz(-0.98139757) q[3];
sx q[3];
rz(-1.792255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.855298) q[0];
sx q[0];
rz(-1.3110302) q[0];
sx q[0];
rz(1.1899765) q[0];
rz(-1.9189934) q[1];
sx q[1];
rz(-1.906955) q[1];
sx q[1];
rz(-2.55866) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8609498) q[0];
sx q[0];
rz(-1.1034729) q[0];
sx q[0];
rz(-1.1790465) q[0];
x q[1];
rz(-3.1276064) q[2];
sx q[2];
rz(-0.99726935) q[2];
sx q[2];
rz(1.7443466) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1057824) q[1];
sx q[1];
rz(-1.3329734) q[1];
sx q[1];
rz(1.4877122) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63151367) q[3];
sx q[3];
rz(-1.9296292) q[3];
sx q[3];
rz(1.7874908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0827834) q[2];
sx q[2];
rz(-1.2502111) q[2];
sx q[2];
rz(1.0397376) q[2];
rz(0.58935634) q[3];
sx q[3];
rz(-1.4732692) q[3];
sx q[3];
rz(1.0065494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37650087) q[0];
sx q[0];
rz(-2.9917175) q[0];
sx q[0];
rz(1.3366706) q[0];
rz(2.916015) q[1];
sx q[1];
rz(-1.0712737) q[1];
sx q[1];
rz(-2.371686) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63084376) q[0];
sx q[0];
rz(-2.0389557) q[0];
sx q[0];
rz(1.7450784) q[0];
x q[1];
rz(0.628148) q[2];
sx q[2];
rz(-0.9315812) q[2];
sx q[2];
rz(-1.01835) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0437255) q[1];
sx q[1];
rz(-1.7683523) q[1];
sx q[1];
rz(1.6476589) q[1];
rz(-3.0287241) q[3];
sx q[3];
rz(-2.5333571) q[3];
sx q[3];
rz(-3.1361767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5169107) q[2];
sx q[2];
rz(-1.3930438) q[2];
sx q[2];
rz(-2.1262271) q[2];
rz(-2.9376302) q[3];
sx q[3];
rz(-1.5927619) q[3];
sx q[3];
rz(-1.7502194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2400951) q[0];
sx q[0];
rz(-0.63145852) q[0];
sx q[0];
rz(-2.748306) q[0];
rz(0.45010629) q[1];
sx q[1];
rz(-1.7186807) q[1];
sx q[1];
rz(0.030287655) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5092376) q[0];
sx q[0];
rz(-1.3233174) q[0];
sx q[0];
rz(-1.6028274) q[0];
rz(-pi) q[1];
rz(2.8856336) q[2];
sx q[2];
rz(-0.82590196) q[2];
sx q[2];
rz(1.2307375) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4612164) q[1];
sx q[1];
rz(-0.54619563) q[1];
sx q[1];
rz(3.097018) q[1];
rz(-pi) q[2];
rz(1.7025331) q[3];
sx q[3];
rz(-1.3088969) q[3];
sx q[3];
rz(-1.9528263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1071189) q[2];
sx q[2];
rz(-1.5849042) q[2];
sx q[2];
rz(-2.0195473) q[2];
rz(-1.0551039) q[3];
sx q[3];
rz(-2.7542346) q[3];
sx q[3];
rz(1.25157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17290641) q[0];
sx q[0];
rz(-2.1538669) q[0];
sx q[0];
rz(-0.77447844) q[0];
rz(-0.6340181) q[1];
sx q[1];
rz(-2.1114383) q[1];
sx q[1];
rz(-1.0672306) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9241079) q[0];
sx q[0];
rz(-1.587968) q[0];
sx q[0];
rz(0.73826684) q[0];
rz(-1.2634694) q[2];
sx q[2];
rz(-1.5241703) q[2];
sx q[2];
rz(-1.4263375) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6065403) q[1];
sx q[1];
rz(-2.7379704) q[1];
sx q[1];
rz(0.89151793) q[1];
rz(-pi) q[2];
rz(-1.0916294) q[3];
sx q[3];
rz(-1.7288343) q[3];
sx q[3];
rz(-0.8534067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.60712236) q[2];
sx q[2];
rz(-3.0202713) q[2];
sx q[2];
rz(-1.8758476) q[2];
rz(0.040146116) q[3];
sx q[3];
rz(-2.7669192) q[3];
sx q[3];
rz(-0.039073959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6604615) q[0];
sx q[0];
rz(-0.88812319) q[0];
sx q[0];
rz(1.2598502) q[0];
rz(-1.4866359) q[1];
sx q[1];
rz(-2.1075552) q[1];
sx q[1];
rz(0.0025509603) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77230763) q[0];
sx q[0];
rz(-2.9678759) q[0];
sx q[0];
rz(0.80887167) q[0];
rz(-1.4223406) q[2];
sx q[2];
rz(-3.0768148) q[2];
sx q[2];
rz(1.5671519) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.59460179) q[1];
sx q[1];
rz(-0.78557116) q[1];
sx q[1];
rz(2.5681096) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6278159) q[3];
sx q[3];
rz(-2.7064763) q[3];
sx q[3];
rz(-0.42729898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.5551787) q[2];
sx q[2];
rz(-1.6224344) q[2];
sx q[2];
rz(-2.4671538) q[2];
rz(-2.0241578) q[3];
sx q[3];
rz(-2.1148465) q[3];
sx q[3];
rz(-2.8286772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085070327) q[0];
sx q[0];
rz(-1.1310348) q[0];
sx q[0];
rz(-0.52666589) q[0];
rz(-0.3479192) q[1];
sx q[1];
rz(-1.2384474) q[1];
sx q[1];
rz(1.7927982) q[1];
rz(0.5402901) q[2];
sx q[2];
rz(-0.39763309) q[2];
sx q[2];
rz(2.8297505) q[2];
rz(-1.717129) q[3];
sx q[3];
rz(-1.5884593) q[3];
sx q[3];
rz(-2.9149169) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
