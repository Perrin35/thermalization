OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4797526) q[0];
sx q[0];
rz(-2.2979484) q[0];
sx q[0];
rz(2.9736829) q[0];
rz(1.1711988) q[1];
sx q[1];
rz(3.436915) q[1];
sx q[1];
rz(9.480939) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1355609) q[0];
sx q[0];
rz(-0.52838415) q[0];
sx q[0];
rz(1.1392659) q[0];
rz(1.2916318) q[2];
sx q[2];
rz(-0.66510519) q[2];
sx q[2];
rz(1.6601738) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9819298) q[1];
sx q[1];
rz(-0.61239457) q[1];
sx q[1];
rz(1.0931404) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5264741) q[3];
sx q[3];
rz(-1.830415) q[3];
sx q[3];
rz(-2.0093902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7636259) q[2];
sx q[2];
rz(-0.28181919) q[2];
sx q[2];
rz(-2.7089233) q[2];
rz(-1.9487322) q[3];
sx q[3];
rz(-1.9038707) q[3];
sx q[3];
rz(-0.38309923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52779657) q[0];
sx q[0];
rz(-0.48848099) q[0];
sx q[0];
rz(1.312785) q[0];
rz(-2.9361172) q[1];
sx q[1];
rz(-0.97646362) q[1];
sx q[1];
rz(1.9899433) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8173556) q[0];
sx q[0];
rz(-1.2756057) q[0];
sx q[0];
rz(2.2833707) q[0];
rz(-2.5580514) q[2];
sx q[2];
rz(-1.984664) q[2];
sx q[2];
rz(-1.5577424) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.68576751) q[1];
sx q[1];
rz(-2.0058504) q[1];
sx q[1];
rz(-2.0352092) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50206708) q[3];
sx q[3];
rz(-1.2289398) q[3];
sx q[3];
rz(-1.348192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1318704) q[2];
sx q[2];
rz(-1.4902318) q[2];
sx q[2];
rz(-2.9197664) q[2];
rz(2.7644073) q[3];
sx q[3];
rz(-2.714034) q[3];
sx q[3];
rz(0.77243531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8310228) q[0];
sx q[0];
rz(-0.092386827) q[0];
sx q[0];
rz(0.036852766) q[0];
rz(0.82551461) q[1];
sx q[1];
rz(-1.8258391) q[1];
sx q[1];
rz(0.056578606) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.409257) q[0];
sx q[0];
rz(-2.0154672) q[0];
sx q[0];
rz(0.94629855) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47358863) q[2];
sx q[2];
rz(-2.8550672) q[2];
sx q[2];
rz(-2.5055656) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.50443447) q[1];
sx q[1];
rz(-1.27099) q[1];
sx q[1];
rz(1.3624886) q[1];
rz(2.644396) q[3];
sx q[3];
rz(-2.433784) q[3];
sx q[3];
rz(-1.4149815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.27292192) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(-2.2154714) q[2];
rz(-0.55666322) q[3];
sx q[3];
rz(-2.848048) q[3];
sx q[3];
rz(2.0986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.91519231) q[0];
sx q[0];
rz(-2.4214348) q[0];
sx q[0];
rz(2.3994989) q[0];
rz(-1.1391976) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(0.46359584) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6620561) q[0];
sx q[0];
rz(-2.2502796) q[0];
sx q[0];
rz(0.38328538) q[0];
rz(-pi) q[1];
rz(-2.8088403) q[2];
sx q[2];
rz(-1.5126265) q[2];
sx q[2];
rz(1.0156877) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5220118) q[1];
sx q[1];
rz(-1.3995692) q[1];
sx q[1];
rz(-0.79230688) q[1];
rz(2.2673244) q[3];
sx q[3];
rz(-2.4618751) q[3];
sx q[3];
rz(0.10248871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3670369) q[2];
sx q[2];
rz(-2.98428) q[2];
sx q[2];
rz(0.049499361) q[2];
rz(-3.0130623) q[3];
sx q[3];
rz(-1.5934207) q[3];
sx q[3];
rz(-3.0226504) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13609919) q[0];
sx q[0];
rz(-0.43731421) q[0];
sx q[0];
rz(0.29770011) q[0];
rz(-2.659335) q[1];
sx q[1];
rz(-2.3869956) q[1];
sx q[1];
rz(-0.94435) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7992295) q[0];
sx q[0];
rz(-2.9276597) q[0];
sx q[0];
rz(-1.3571204) q[0];
rz(-0.77563939) q[2];
sx q[2];
rz(-1.8619985) q[2];
sx q[2];
rz(1.8679384) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.22630616) q[1];
sx q[1];
rz(-1.5686791) q[1];
sx q[1];
rz(1.5096942) q[1];
x q[2];
rz(-2.6258351) q[3];
sx q[3];
rz(-0.48832794) q[3];
sx q[3];
rz(0.26667903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.258761) q[2];
sx q[2];
rz(-1.0927039) q[2];
sx q[2];
rz(3.0333701) q[2];
rz(-0.0023068874) q[3];
sx q[3];
rz(-1.6121515) q[3];
sx q[3];
rz(2.8172857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7047983) q[0];
sx q[0];
rz(-2.695485) q[0];
sx q[0];
rz(-2.5571402) q[0];
rz(-2.2553518) q[1];
sx q[1];
rz(-2.5247572) q[1];
sx q[1];
rz(-3.086673) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6829837) q[0];
sx q[0];
rz(-2.8603472) q[0];
sx q[0];
rz(-1.8722697) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5149649) q[2];
sx q[2];
rz(-2.8179114) q[2];
sx q[2];
rz(1.8813546) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.9565441) q[1];
sx q[1];
rz(-1.1951606) q[1];
sx q[1];
rz(-2.5513785) q[1];
rz(-pi) q[2];
rz(-2.6782126) q[3];
sx q[3];
rz(-1.0372835) q[3];
sx q[3];
rz(-0.86138553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3871258) q[2];
sx q[2];
rz(-0.12067623) q[2];
sx q[2];
rz(2.1248655) q[2];
rz(-2.5975442) q[3];
sx q[3];
rz(-2.7816911) q[3];
sx q[3];
rz(1.3325161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.465437) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(-0.28453919) q[0];
rz(0.94447213) q[1];
sx q[1];
rz(-1.9453134) q[1];
sx q[1];
rz(-2.231266) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0383354) q[0];
sx q[0];
rz(-0.084780134) q[0];
sx q[0];
rz(-2.2708562) q[0];
rz(-pi) q[1];
rz(-1.9867284) q[2];
sx q[2];
rz(-1.8850733) q[2];
sx q[2];
rz(1.8679801) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.42156223) q[1];
sx q[1];
rz(-1.3898464) q[1];
sx q[1];
rz(0.49444316) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1479285) q[3];
sx q[3];
rz(-2.2054407) q[3];
sx q[3];
rz(1.5594547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3900782) q[2];
sx q[2];
rz(-3.0780767) q[2];
sx q[2];
rz(2.2195623) q[2];
rz(2.5743124) q[3];
sx q[3];
rz(-1.7276948) q[3];
sx q[3];
rz(2.1218307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2475125) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(3.0122053) q[0];
rz(-0.63240504) q[1];
sx q[1];
rz(-1.0267195) q[1];
sx q[1];
rz(-2.8410889) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6209517) q[0];
sx q[0];
rz(-1.2004939) q[0];
sx q[0];
rz(-2.305549) q[0];
rz(2.1883165) q[2];
sx q[2];
rz(-2.2308908) q[2];
sx q[2];
rz(-2.2613139) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6943372) q[1];
sx q[1];
rz(-1.9112504) q[1];
sx q[1];
rz(-1.7561595) q[1];
x q[2];
rz(-2.3281906) q[3];
sx q[3];
rz(-1.8040931) q[3];
sx q[3];
rz(2.5113311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.58632103) q[2];
sx q[2];
rz(-2.1466612) q[2];
sx q[2];
rz(2.3596181) q[2];
rz(0.55109763) q[3];
sx q[3];
rz(-1.7588153) q[3];
sx q[3];
rz(-0.15792318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5732116) q[0];
sx q[0];
rz(-1.9848354) q[0];
sx q[0];
rz(-3.0138299) q[0];
rz(2.5993775) q[1];
sx q[1];
rz(-0.95710373) q[1];
sx q[1];
rz(-2.382747) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.471506) q[0];
sx q[0];
rz(-0.42450464) q[0];
sx q[0];
rz(1.5295117) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60894572) q[2];
sx q[2];
rz(-1.7296089) q[2];
sx q[2];
rz(1.2283404) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2821741) q[1];
sx q[1];
rz(-0.71600435) q[1];
sx q[1];
rz(0.23846682) q[1];
x q[2];
rz(-1.5522478) q[3];
sx q[3];
rz(-2.4759001) q[3];
sx q[3];
rz(1.1902283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0163429) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(-2.6515567) q[2];
rz(1.4222493) q[3];
sx q[3];
rz(-1.9336721) q[3];
sx q[3];
rz(-2.0786044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35995099) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(-0.075335659) q[0];
rz(-2.244859) q[1];
sx q[1];
rz(-1.9672085) q[1];
sx q[1];
rz(0.60992253) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88203726) q[0];
sx q[0];
rz(-1.2872739) q[0];
sx q[0];
rz(2.3638704) q[0];
rz(1.8462734) q[2];
sx q[2];
rz(-1.9343978) q[2];
sx q[2];
rz(1.3678577) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6590609) q[1];
sx q[1];
rz(-2.3606803) q[1];
sx q[1];
rz(-0.34300967) q[1];
rz(-1.0697332) q[3];
sx q[3];
rz(-0.89322972) q[3];
sx q[3];
rz(-0.60615218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23218368) q[2];
sx q[2];
rz(-2.3164618) q[2];
sx q[2];
rz(-0.71371901) q[2];
rz(2.7632726) q[3];
sx q[3];
rz(-2.648073) q[3];
sx q[3];
rz(2.2617214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.338035) q[0];
sx q[0];
rz(-1.9914347) q[0];
sx q[0];
rz(1.5557355) q[0];
rz(-0.65080416) q[1];
sx q[1];
rz(-1.4918068) q[1];
sx q[1];
rz(3.0202958) q[1];
rz(-1.3748319) q[2];
sx q[2];
rz(-1.5407731) q[2];
sx q[2];
rz(-2.4827448) q[2];
rz(2.0416904) q[3];
sx q[3];
rz(-1.7970016) q[3];
sx q[3];
rz(-2.5947528) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
