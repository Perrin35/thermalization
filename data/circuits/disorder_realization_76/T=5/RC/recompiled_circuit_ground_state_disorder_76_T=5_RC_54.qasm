OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6948833) q[0];
sx q[0];
rz(-1.0816242) q[0];
sx q[0];
rz(-0.73944902) q[0];
rz(-0.75955716) q[1];
sx q[1];
rz(-1.8269202) q[1];
sx q[1];
rz(-1.4256328) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6524871) q[0];
sx q[0];
rz(-1.28363) q[0];
sx q[0];
rz(2.1225568) q[0];
rz(-pi) q[1];
rz(1.8427467) q[2];
sx q[2];
rz(-1.0894948) q[2];
sx q[2];
rz(1.1084194) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.718218) q[1];
sx q[1];
rz(-1.4932695) q[1];
sx q[1];
rz(-1.422387) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0979068) q[3];
sx q[3];
rz(-2.2891392) q[3];
sx q[3];
rz(-1.1506611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2322959) q[2];
sx q[2];
rz(-2.2984419) q[2];
sx q[2];
rz(-0.24093957) q[2];
rz(3.1124034) q[3];
sx q[3];
rz(-1.3386644) q[3];
sx q[3];
rz(1.1791112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1643243) q[0];
sx q[0];
rz(-1.203953) q[0];
sx q[0];
rz(-0.48467317) q[0];
rz(2.4618705) q[1];
sx q[1];
rz(-1.2815963) q[1];
sx q[1];
rz(-1.9814804) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.401448) q[0];
sx q[0];
rz(-0.30456802) q[0];
sx q[0];
rz(-3.079971) q[0];
rz(-pi) q[1];
rz(-0.76267879) q[2];
sx q[2];
rz(-0.87717036) q[2];
sx q[2];
rz(0.1887624) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1801123) q[1];
sx q[1];
rz(-0.65458502) q[1];
sx q[1];
rz(-2.727319) q[1];
rz(-0.53547041) q[3];
sx q[3];
rz(-2.5512085) q[3];
sx q[3];
rz(-1.2108506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.29740563) q[2];
sx q[2];
rz(-0.30218267) q[2];
sx q[2];
rz(-0.45787946) q[2];
rz(2.0186021) q[3];
sx q[3];
rz(-0.99959683) q[3];
sx q[3];
rz(2.5751953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26894012) q[0];
sx q[0];
rz(-2.432423) q[0];
sx q[0];
rz(-3.1100682) q[0];
rz(-2.8541376) q[1];
sx q[1];
rz(-0.87702409) q[1];
sx q[1];
rz(-1.215975) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4327964) q[0];
sx q[0];
rz(-0.93892083) q[0];
sx q[0];
rz(0.89544501) q[0];
rz(-pi) q[1];
rz(-0.48451938) q[2];
sx q[2];
rz(-1.7555825) q[2];
sx q[2];
rz(-2.4444524) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.79920125) q[1];
sx q[1];
rz(-1.2974129) q[1];
sx q[1];
rz(-0.28634109) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4803355) q[3];
sx q[3];
rz(-1.0399858) q[3];
sx q[3];
rz(-2.9531933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.03881255) q[2];
sx q[2];
rz(-1.5993885) q[2];
sx q[2];
rz(2.3056324) q[2];
rz(1.7838259) q[3];
sx q[3];
rz(-0.72503763) q[3];
sx q[3];
rz(-0.74762216) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2826071) q[0];
sx q[0];
rz(-1.7361807) q[0];
sx q[0];
rz(-0.080168515) q[0];
rz(1.0824341) q[1];
sx q[1];
rz(-0.23324649) q[1];
sx q[1];
rz(-1.8128043) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4565312) q[0];
sx q[0];
rz(-1.0624806) q[0];
sx q[0];
rz(1.9446745) q[0];
x q[1];
rz(2.6615449) q[2];
sx q[2];
rz(-2.7809395) q[2];
sx q[2];
rz(0.03483054) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3995953) q[1];
sx q[1];
rz(-0.4542225) q[1];
sx q[1];
rz(1.9059559) q[1];
rz(-pi) q[2];
rz(-0.53502632) q[3];
sx q[3];
rz(-0.38540977) q[3];
sx q[3];
rz(-0.42675323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7666011) q[2];
sx q[2];
rz(-0.98321715) q[2];
sx q[2];
rz(0.90744606) q[2];
rz(2.4713016) q[3];
sx q[3];
rz(-2.3136316) q[3];
sx q[3];
rz(-1.7214187) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2403253) q[0];
sx q[0];
rz(-2.2690052) q[0];
sx q[0];
rz(-2.7101044) q[0];
rz(-2.0101428) q[1];
sx q[1];
rz(-1.1791469) q[1];
sx q[1];
rz(-0.44050899) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9517072) q[0];
sx q[0];
rz(-2.0312956) q[0];
sx q[0];
rz(-3.0104464) q[0];
rz(-pi) q[1];
rz(-1.9856057) q[2];
sx q[2];
rz(-1.3813263) q[2];
sx q[2];
rz(-2.1304325) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.1738834) q[1];
sx q[1];
rz(-1.6102127) q[1];
sx q[1];
rz(-0.33052488) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93705658) q[3];
sx q[3];
rz(-2.0645294) q[3];
sx q[3];
rz(2.1449094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.12209192) q[2];
sx q[2];
rz(-0.92635218) q[2];
sx q[2];
rz(2.7276373) q[2];
rz(1.3881989) q[3];
sx q[3];
rz(-1.5940758) q[3];
sx q[3];
rz(-3.0164914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6269161) q[0];
sx q[0];
rz(-1.7358945) q[0];
sx q[0];
rz(2.2077014) q[0];
rz(-2.4712708) q[1];
sx q[1];
rz(-1.2443845) q[1];
sx q[1];
rz(-1.3353039) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2349699) q[0];
sx q[0];
rz(-2.1701394) q[0];
sx q[0];
rz(-0.12110981) q[0];
rz(-pi) q[1];
rz(2.983903) q[2];
sx q[2];
rz(-1.5458428) q[2];
sx q[2];
rz(-1.4098997) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2833497) q[1];
sx q[1];
rz(-1.5482386) q[1];
sx q[1];
rz(1.683504) q[1];
rz(0.63702668) q[3];
sx q[3];
rz(-1.5638132) q[3];
sx q[3];
rz(-1.03656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.89207092) q[2];
sx q[2];
rz(-2.8688909) q[2];
sx q[2];
rz(-1.8708694) q[2];
rz(-1.7719841) q[3];
sx q[3];
rz(-1.2415875) q[3];
sx q[3];
rz(2.6197267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9524566) q[0];
sx q[0];
rz(-0.58681762) q[0];
sx q[0];
rz(2.9879046) q[0];
rz(2.9391089) q[1];
sx q[1];
rz(-1.9053562) q[1];
sx q[1];
rz(-0.94857803) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0517387) q[0];
sx q[0];
rz(-1.7719381) q[0];
sx q[0];
rz(0.62369831) q[0];
rz(-pi) q[1];
rz(-2.8178627) q[2];
sx q[2];
rz(-0.95697953) q[2];
sx q[2];
rz(-2.987189) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3707917) q[1];
sx q[1];
rz(-1.3389412) q[1];
sx q[1];
rz(-1.8370017) q[1];
x q[2];
rz(-0.13410577) q[3];
sx q[3];
rz(-0.92536345) q[3];
sx q[3];
rz(-2.9112618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1367246) q[2];
sx q[2];
rz(-1.0978881) q[2];
sx q[2];
rz(-0.0014121545) q[2];
rz(-3.0794365) q[3];
sx q[3];
rz(-1.4096189) q[3];
sx q[3];
rz(-2.1803161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75981265) q[0];
sx q[0];
rz(-2.3842922) q[0];
sx q[0];
rz(1.9267474) q[0];
rz(-0.75792056) q[1];
sx q[1];
rz(-0.52752033) q[1];
sx q[1];
rz(-0.55799276) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3787631) q[0];
sx q[0];
rz(-0.42233322) q[0];
sx q[0];
rz(2.2076553) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5900389) q[2];
sx q[2];
rz(-1.6948912) q[2];
sx q[2];
rz(0.53260224) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.87822014) q[1];
sx q[1];
rz(-1.0910631) q[1];
sx q[1];
rz(0.31444506) q[1];
rz(-pi) q[2];
rz(-1.8479061) q[3];
sx q[3];
rz(-1.596611) q[3];
sx q[3];
rz(1.1518948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5743635) q[2];
sx q[2];
rz(-1.1939253) q[2];
sx q[2];
rz(-2.8723259) q[2];
rz(-1.1632129) q[3];
sx q[3];
rz(-0.60267699) q[3];
sx q[3];
rz(0.24063024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6729386) q[0];
sx q[0];
rz(-2.1042295) q[0];
sx q[0];
rz(-1.0733676) q[0];
rz(2.0138373) q[1];
sx q[1];
rz(-2.914371) q[1];
sx q[1];
rz(-2.5849297) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23281413) q[0];
sx q[0];
rz(-1.4289843) q[0];
sx q[0];
rz(-1.892426) q[0];
rz(-0.78293856) q[2];
sx q[2];
rz(-0.73854827) q[2];
sx q[2];
rz(-2.1434181) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8175259) q[1];
sx q[1];
rz(-2.3285667) q[1];
sx q[1];
rz(2.6387385) q[1];
x q[2];
rz(0.23454097) q[3];
sx q[3];
rz(-1.4734771) q[3];
sx q[3];
rz(-0.095712599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.90159455) q[2];
sx q[2];
rz(-1.0866714) q[2];
sx q[2];
rz(-1.979801) q[2];
rz(-2.6084172) q[3];
sx q[3];
rz(-0.52459255) q[3];
sx q[3];
rz(0.1575135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6530782) q[0];
sx q[0];
rz(-2.1670659) q[0];
sx q[0];
rz(0.59481204) q[0];
rz(-0.57890233) q[1];
sx q[1];
rz(-1.4756823) q[1];
sx q[1];
rz(-0.65931177) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66856495) q[0];
sx q[0];
rz(-0.94946948) q[0];
sx q[0];
rz(-1.4061808) q[0];
rz(-pi) q[1];
rz(1.1289146) q[2];
sx q[2];
rz(-2.127248) q[2];
sx q[2];
rz(2.93612) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.51344556) q[1];
sx q[1];
rz(-1.7127348) q[1];
sx q[1];
rz(-1.3238841) q[1];
rz(-pi) q[2];
rz(2.8774977) q[3];
sx q[3];
rz(-1.8349097) q[3];
sx q[3];
rz(-2.7434289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.79002964) q[2];
sx q[2];
rz(-1.1976676) q[2];
sx q[2];
rz(-3.0885546) q[2];
rz(1.7985581) q[3];
sx q[3];
rz(-1.8648632) q[3];
sx q[3];
rz(-3.067335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85528436) q[0];
sx q[0];
rz(-1.3636148) q[0];
sx q[0];
rz(-1.6028945) q[0];
rz(0.098943624) q[1];
sx q[1];
rz(-1.3700486) q[1];
sx q[1];
rz(0.055421967) q[1];
rz(1.5621875) q[2];
sx q[2];
rz(-1.1653524) q[2];
sx q[2];
rz(-0.63962519) q[2];
rz(1.4138447) q[3];
sx q[3];
rz(-2.3060287) q[3];
sx q[3];
rz(-2.9097478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
