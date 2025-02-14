OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8813397) q[0];
sx q[0];
rz(-0.94085675) q[0];
sx q[0];
rz(2.9139304) q[0];
rz(-2.8582299) q[1];
sx q[1];
rz(-0.41937399) q[1];
sx q[1];
rz(1.286932) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0193699) q[0];
sx q[0];
rz(-1.625018) q[0];
sx q[0];
rz(-2.0383459) q[0];
rz(2.020316) q[2];
sx q[2];
rz(-1.9971022) q[2];
sx q[2];
rz(0.21838494) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.097295105) q[1];
sx q[1];
rz(-1.0781647) q[1];
sx q[1];
rz(0.45483847) q[1];
rz(-1.257913) q[3];
sx q[3];
rz(-2.1863424) q[3];
sx q[3];
rz(2.2039977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8925573) q[2];
sx q[2];
rz(-1.0113357) q[2];
sx q[2];
rz(-2.2300143) q[2];
rz(-0.75561953) q[3];
sx q[3];
rz(-0.32114649) q[3];
sx q[3];
rz(0.49629456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3590473) q[0];
sx q[0];
rz(-1.8308715) q[0];
sx q[0];
rz(-2.0027335) q[0];
rz(0.62659872) q[1];
sx q[1];
rz(-2.7732924) q[1];
sx q[1];
rz(-1.0777333) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8846466) q[0];
sx q[0];
rz(-2.4532689) q[0];
sx q[0];
rz(-0.55082816) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1280649) q[2];
sx q[2];
rz(-1.2473543) q[2];
sx q[2];
rz(-1.5087104) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5087532) q[1];
sx q[1];
rz(-1.5037854) q[1];
sx q[1];
rz(-2.8772023) q[1];
rz(-pi) q[2];
rz(2.2935981) q[3];
sx q[3];
rz(-1.2395596) q[3];
sx q[3];
rz(-0.02324638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.44800147) q[2];
sx q[2];
rz(-2.3841264) q[2];
sx q[2];
rz(2.9288911) q[2];
rz(1.5564144) q[3];
sx q[3];
rz(-2.3978265) q[3];
sx q[3];
rz(-1.0630382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0001275) q[0];
sx q[0];
rz(-0.26697049) q[0];
sx q[0];
rz(-2.2947327) q[0];
rz(2.8894539) q[1];
sx q[1];
rz(-2.3138901) q[1];
sx q[1];
rz(-0.65021461) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92714053) q[0];
sx q[0];
rz(-0.34159476) q[0];
sx q[0];
rz(0.13169698) q[0];
rz(-pi) q[1];
rz(-2.3450801) q[2];
sx q[2];
rz(-1.0625417) q[2];
sx q[2];
rz(-0.69895335) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0421335) q[1];
sx q[1];
rz(-0.317527) q[1];
sx q[1];
rz(-1.9327764) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3350983) q[3];
sx q[3];
rz(-1.8915081) q[3];
sx q[3];
rz(2.210828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.83027679) q[2];
sx q[2];
rz(-0.69861424) q[2];
sx q[2];
rz(-0.96381956) q[2];
rz(1.7802995) q[3];
sx q[3];
rz(-2.6908974) q[3];
sx q[3];
rz(3.0986339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2394102) q[0];
sx q[0];
rz(-0.22265156) q[0];
sx q[0];
rz(-1.0182925) q[0];
rz(0.88515431) q[1];
sx q[1];
rz(-2.8926909) q[1];
sx q[1];
rz(0.28908602) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24482306) q[0];
sx q[0];
rz(-2.6894248) q[0];
sx q[0];
rz(-0.74613476) q[0];
x q[1];
rz(-2.5968005) q[2];
sx q[2];
rz(-2.4281686) q[2];
sx q[2];
rz(-0.95470828) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8377467) q[1];
sx q[1];
rz(-2.0833384) q[1];
sx q[1];
rz(-0.20767494) q[1];
rz(-pi) q[2];
rz(-2.9226275) q[3];
sx q[3];
rz(-1.4462399) q[3];
sx q[3];
rz(-1.3606461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3503512) q[2];
sx q[2];
rz(-1.6356607) q[2];
sx q[2];
rz(-1.2887456) q[2];
rz(-0.57981235) q[3];
sx q[3];
rz(-0.47477397) q[3];
sx q[3];
rz(-0.60540664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30699214) q[0];
sx q[0];
rz(-2.6267316) q[0];
sx q[0];
rz(-0.32421625) q[0];
rz(-0.038837198) q[1];
sx q[1];
rz(-0.7380929) q[1];
sx q[1];
rz(2.0447581) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0457458) q[0];
sx q[0];
rz(-0.37558324) q[0];
sx q[0];
rz(2.4235241) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4154424) q[2];
sx q[2];
rz(-2.5553779) q[2];
sx q[2];
rz(1.6833351) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3009764) q[1];
sx q[1];
rz(-0.30472091) q[1];
sx q[1];
rz(-1.4304377) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.092336307) q[3];
sx q[3];
rz(-1.9895377) q[3];
sx q[3];
rz(3.0224722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0754806) q[2];
sx q[2];
rz(-2.8754063) q[2];
sx q[2];
rz(0.038662635) q[2];
rz(1.4085116) q[3];
sx q[3];
rz(-1.446529) q[3];
sx q[3];
rz(-2.8601638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.6283145) q[0];
sx q[0];
rz(-2.3230041) q[0];
sx q[0];
rz(-2.1224838) q[0];
rz(0.45711532) q[1];
sx q[1];
rz(-0.26908427) q[1];
sx q[1];
rz(1.4310744) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5810878) q[0];
sx q[0];
rz(-1.3700587) q[0];
sx q[0];
rz(-1.5439347) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.25361714) q[2];
sx q[2];
rz(-2.1448958) q[2];
sx q[2];
rz(-0.92437896) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0455605) q[1];
sx q[1];
rz(-1.9304196) q[1];
sx q[1];
rz(-0.093861967) q[1];
rz(-2.01675) q[3];
sx q[3];
rz(-0.77174458) q[3];
sx q[3];
rz(1.2260557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9223601) q[2];
sx q[2];
rz(-1.7722426) q[2];
sx q[2];
rz(0.57979453) q[2];
rz(0.29948768) q[3];
sx q[3];
rz(-0.53286415) q[3];
sx q[3];
rz(-1.2286435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2806468) q[0];
sx q[0];
rz(-1.4806643) q[0];
sx q[0];
rz(-3.0378367) q[0];
rz(-0.62458986) q[1];
sx q[1];
rz(-2.2207405) q[1];
sx q[1];
rz(-1.7594899) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64829761) q[0];
sx q[0];
rz(-1.7837423) q[0];
sx q[0];
rz(0.84342028) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8417009) q[2];
sx q[2];
rz(-1.4647398) q[2];
sx q[2];
rz(2.7404355) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.15006062) q[1];
sx q[1];
rz(-2.074291) q[1];
sx q[1];
rz(0.30725422) q[1];
x q[2];
rz(0.35993536) q[3];
sx q[3];
rz(-2.4483134) q[3];
sx q[3];
rz(0.14132796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.78911191) q[2];
sx q[2];
rz(-1.666297) q[2];
sx q[2];
rz(2.2769807) q[2];
rz(0.23884808) q[3];
sx q[3];
rz(-2.3819203) q[3];
sx q[3];
rz(0.71340942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(0.1777986) q[0];
sx q[0];
rz(-0.56121427) q[0];
sx q[0];
rz(0.23502769) q[0];
rz(0.66559732) q[1];
sx q[1];
rz(-2.0033629) q[1];
sx q[1];
rz(1.4733018) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.145133) q[0];
sx q[0];
rz(-1.0211103) q[0];
sx q[0];
rz(-2.8966146) q[0];
rz(-pi) q[1];
rz(-2.8741415) q[2];
sx q[2];
rz(-2.6595734) q[2];
sx q[2];
rz(-2.6798525) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7104147) q[1];
sx q[1];
rz(-1.7542375) q[1];
sx q[1];
rz(1.4967493) q[1];
x q[2];
rz(-3.0489122) q[3];
sx q[3];
rz(-0.4527992) q[3];
sx q[3];
rz(-1.6860675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9792446) q[2];
sx q[2];
rz(-0.20077106) q[2];
sx q[2];
rz(2.6806504) q[2];
rz(-1.0848684) q[3];
sx q[3];
rz(-2.3960787) q[3];
sx q[3];
rz(2.357024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0775065) q[0];
sx q[0];
rz(-0.53090799) q[0];
sx q[0];
rz(0.60894668) q[0];
rz(-2.5755836) q[1];
sx q[1];
rz(-1.1606263) q[1];
sx q[1];
rz(-1.2014679) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5574359) q[0];
sx q[0];
rz(-3.0058001) q[0];
sx q[0];
rz(1.6312172) q[0];
rz(-pi) q[1];
rz(0.20862357) q[2];
sx q[2];
rz(-1.0917328) q[2];
sx q[2];
rz(-0.63392679) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4200099) q[1];
sx q[1];
rz(-2.3370565) q[1];
sx q[1];
rz(-1.700042) q[1];
rz(-pi) q[2];
rz(1.8747599) q[3];
sx q[3];
rz(-2.0018775) q[3];
sx q[3];
rz(1.4747185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0386049) q[2];
sx q[2];
rz(-1.0706341) q[2];
sx q[2];
rz(-1.1919682) q[2];
rz(-1.0209171) q[3];
sx q[3];
rz(-2.9396785) q[3];
sx q[3];
rz(-1.6010223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22333071) q[0];
sx q[0];
rz(-2.9409565) q[0];
sx q[0];
rz(2.1740792) q[0];
rz(1.7009023) q[1];
sx q[1];
rz(-0.43165019) q[1];
sx q[1];
rz(-0.38223019) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4995196) q[0];
sx q[0];
rz(-1.0040305) q[0];
sx q[0];
rz(2.9261409) q[0];
x q[1];
rz(1.9873496) q[2];
sx q[2];
rz(-1.9723411) q[2];
sx q[2];
rz(2.7018715) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9542071) q[1];
sx q[1];
rz(-1.6897196) q[1];
sx q[1];
rz(-2.3559761) q[1];
x q[2];
rz(-0.02703826) q[3];
sx q[3];
rz(-0.54618764) q[3];
sx q[3];
rz(2.9423713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.050015673) q[2];
sx q[2];
rz(-0.86805582) q[2];
sx q[2];
rz(2.3559605) q[2];
rz(0.026963726) q[3];
sx q[3];
rz(-1.5188768) q[3];
sx q[3];
rz(-0.54076076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80326573) q[0];
sx q[0];
rz(-1.4563518) q[0];
sx q[0];
rz(-0.8932054) q[0];
rz(0.66435736) q[1];
sx q[1];
rz(-1.9021481) q[1];
sx q[1];
rz(-0.93217168) q[1];
rz(2.4128466) q[2];
sx q[2];
rz(-1.797429) q[2];
sx q[2];
rz(0.28354473) q[2];
rz(3.1278004) q[3];
sx q[3];
rz(-2.5368555) q[3];
sx q[3];
rz(0.93750877) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
