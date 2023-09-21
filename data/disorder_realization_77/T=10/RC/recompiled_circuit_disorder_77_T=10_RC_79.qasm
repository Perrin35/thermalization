OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(-0.84364426) q[0];
sx q[0];
rz(0.16790976) q[0];
rz(1.1711988) q[1];
sx q[1];
rz(3.436915) q[1];
sx q[1];
rz(9.480939) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1355609) q[0];
sx q[0];
rz(-0.52838415) q[0];
sx q[0];
rz(2.0023268) q[0];
rz(-2.9287455) q[2];
sx q[2];
rz(-0.93570645) q[2];
sx q[2];
rz(-1.1320621) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4178884) q[1];
sx q[1];
rz(-1.0350409) q[1];
sx q[1];
rz(0.31236155) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25986259) q[3];
sx q[3];
rz(-1.5279603) q[3];
sx q[3];
rz(2.6916137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7636259) q[2];
sx q[2];
rz(-2.8597735) q[2];
sx q[2];
rz(2.7089233) q[2];
rz(-1.9487322) q[3];
sx q[3];
rz(-1.9038707) q[3];
sx q[3];
rz(2.7584934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6137961) q[0];
sx q[0];
rz(-2.6531117) q[0];
sx q[0];
rz(-1.312785) q[0];
rz(0.20547543) q[1];
sx q[1];
rz(-2.165129) q[1];
sx q[1];
rz(1.1516494) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8173556) q[0];
sx q[0];
rz(-1.2756057) q[0];
sx q[0];
rz(-0.8582219) q[0];
rz(-pi) q[1];
rz(-0.67302455) q[2];
sx q[2];
rz(-2.4403799) q[2];
sx q[2];
rz(-2.6075624) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.68576751) q[1];
sx q[1];
rz(-2.0058504) q[1];
sx q[1];
rz(-2.0352092) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1851951) q[3];
sx q[3];
rz(-2.0413627) q[3];
sx q[3];
rz(3.1009931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.0097222086) q[2];
sx q[2];
rz(-1.4902318) q[2];
sx q[2];
rz(-2.9197664) q[2];
rz(0.37718537) q[3];
sx q[3];
rz(-0.42755869) q[3];
sx q[3];
rz(0.77243531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31056988) q[0];
sx q[0];
rz(-3.0492058) q[0];
sx q[0];
rz(3.1047399) q[0];
rz(-2.316078) q[1];
sx q[1];
rz(-1.8258391) q[1];
sx q[1];
rz(0.056578606) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73233561) q[0];
sx q[0];
rz(-1.1261254) q[0];
sx q[0];
rz(-2.1952941) q[0];
rz(-0.25643202) q[2];
sx q[2];
rz(-1.4415381) q[2];
sx q[2];
rz(-1.749922) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.024836) q[1];
sx q[1];
rz(-0.36326888) q[1];
sx q[1];
rz(0.58961745) q[1];
x q[2];
rz(-0.49719663) q[3];
sx q[3];
rz(-2.433784) q[3];
sx q[3];
rz(-1.4149815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.27292192) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(-2.2154714) q[2];
rz(2.5849294) q[3];
sx q[3];
rz(-0.29354468) q[3];
sx q[3];
rz(1.042897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-2.2264003) q[0];
sx q[0];
rz(-0.72015786) q[0];
sx q[0];
rz(-0.74209374) q[0];
rz(2.0023951) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(-2.6779968) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2984943) q[0];
sx q[0];
rz(-1.2756186) q[0];
sx q[0];
rz(2.28736) q[0];
x q[1];
rz(-0.33275231) q[2];
sx q[2];
rz(-1.6289662) q[2];
sx q[2];
rz(-2.1259049) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1222893) q[1];
sx q[1];
rz(-0.79320723) q[1];
sx q[1];
rz(1.3293468) q[1];
rz(-2.1257524) q[3];
sx q[3];
rz(-1.98588) q[3];
sx q[3];
rz(-2.2500028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7745557) q[2];
sx q[2];
rz(-2.98428) q[2];
sx q[2];
rz(-0.049499361) q[2];
rz(-0.1285304) q[3];
sx q[3];
rz(-1.5481719) q[3];
sx q[3];
rz(-3.0226504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0054935) q[0];
sx q[0];
rz(-2.7042784) q[0];
sx q[0];
rz(-2.8438925) q[0];
rz(2.659335) q[1];
sx q[1];
rz(-0.75459701) q[1];
sx q[1];
rz(-0.94435) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56086841) q[0];
sx q[0];
rz(-1.7797911) q[0];
sx q[0];
rz(3.0955549) q[0];
x q[1];
rz(-1.1733426) q[2];
sx q[2];
rz(-2.3059418) q[2];
sx q[2];
rz(-3.1189001) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3443606) q[1];
sx q[1];
rz(-1.6318983) q[1];
sx q[1];
rz(-3.1394715) q[1];
rz(-pi) q[2];
rz(1.827042) q[3];
sx q[3];
rz(-1.1503997) q[3];
sx q[3];
rz(0.30403578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.258761) q[2];
sx q[2];
rz(-2.0488887) q[2];
sx q[2];
rz(0.10822254) q[2];
rz(0.0023068874) q[3];
sx q[3];
rz(-1.6121515) q[3];
sx q[3];
rz(-2.8172857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.7047983) q[0];
sx q[0];
rz(-0.44610766) q[0];
sx q[0];
rz(-0.58445245) q[0];
rz(0.8862409) q[1];
sx q[1];
rz(-2.5247572) q[1];
sx q[1];
rz(0.054919682) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14558218) q[0];
sx q[0];
rz(-1.8390363) q[0];
sx q[0];
rz(-0.085573816) q[0];
rz(-pi) q[1];
rz(-1.6266277) q[2];
sx q[2];
rz(-0.32368127) q[2];
sx q[2];
rz(1.8813546) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7683148) q[1];
sx q[1];
rz(-2.1149966) q[1];
sx q[1];
rz(-2.0139704) q[1];
x q[2];
rz(2.218607) q[3];
sx q[3];
rz(-2.4499948) q[3];
sx q[3];
rz(-1.637961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3871258) q[2];
sx q[2];
rz(-3.0209164) q[2];
sx q[2];
rz(2.1248655) q[2];
rz(0.54404849) q[3];
sx q[3];
rz(-0.35990158) q[3];
sx q[3];
rz(1.8090766) q[3];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6761557) q[0];
sx q[0];
rz(-0.98452079) q[0];
sx q[0];
rz(-2.8570535) q[0];
rz(0.94447213) q[1];
sx q[1];
rz(-1.9453134) q[1];
sx q[1];
rz(-2.231266) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1032573) q[0];
sx q[0];
rz(-3.0568125) q[0];
sx q[0];
rz(-2.2708562) q[0];
rz(1.9867284) q[2];
sx q[2];
rz(-1.2565194) q[2];
sx q[2];
rz(-1.2736125) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.42156223) q[1];
sx q[1];
rz(-1.3898464) q[1];
sx q[1];
rz(0.49444316) q[1];
rz(-pi) q[2];
rz(-0.99366412) q[3];
sx q[3];
rz(-2.2054407) q[3];
sx q[3];
rz(-1.5594547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3900782) q[2];
sx q[2];
rz(-3.0780767) q[2];
sx q[2];
rz(2.2195623) q[2];
rz(0.56728029) q[3];
sx q[3];
rz(-1.4138979) q[3];
sx q[3];
rz(2.1218307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89408016) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(-0.12938736) q[0];
rz(2.5091876) q[1];
sx q[1];
rz(-1.0267195) q[1];
sx q[1];
rz(0.30050373) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7104615) q[0];
sx q[0];
rz(-2.3346402) q[0];
sx q[0];
rz(2.0956844) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6408765) q[2];
sx q[2];
rz(-0.87101988) q[2];
sx q[2];
rz(-1.7388294) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.064425163) q[1];
sx q[1];
rz(-0.38591138) q[1];
sx q[1];
rz(-0.47972958) q[1];
rz(2.3281906) q[3];
sx q[3];
rz(-1.8040931) q[3];
sx q[3];
rz(0.63026159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.58632103) q[2];
sx q[2];
rz(-2.1466612) q[2];
sx q[2];
rz(-2.3596181) q[2];
rz(-0.55109763) q[3];
sx q[3];
rz(-1.7588153) q[3];
sx q[3];
rz(-2.9836695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5732116) q[0];
sx q[0];
rz(-1.9848354) q[0];
sx q[0];
rz(3.0138299) q[0];
rz(-0.54221517) q[1];
sx q[1];
rz(-0.95710373) q[1];
sx q[1];
rz(-2.382747) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7153873) q[0];
sx q[0];
rz(-1.1466768) q[0];
sx q[0];
rz(-0.018652648) q[0];
x q[1];
rz(0.60894572) q[2];
sx q[2];
rz(-1.4119838) q[2];
sx q[2];
rz(1.2283404) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2488333) q[1];
sx q[1];
rz(-1.4151238) q[1];
sx q[1];
rz(-2.4397736) q[1];
x q[2];
rz(-2.2364053) q[3];
sx q[3];
rz(-1.5593411) q[3];
sx q[3];
rz(-0.36597914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1252497) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(0.49003595) q[2];
rz(1.4222493) q[3];
sx q[3];
rz(-1.2079206) q[3];
sx q[3];
rz(2.0786044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7816417) q[0];
sx q[0];
rz(-0.61976969) q[0];
sx q[0];
rz(-3.066257) q[0];
rz(-2.244859) q[1];
sx q[1];
rz(-1.9672085) q[1];
sx q[1];
rz(0.60992253) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7216435) q[0];
sx q[0];
rz(-0.83166612) q[0];
sx q[0];
rz(1.182611) q[0];
x q[1];
rz(-0.37656017) q[2];
sx q[2];
rz(-1.3137523) q[2];
sx q[2];
rz(0.10274796) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.48253179) q[1];
sx q[1];
rz(-2.3606803) q[1];
sx q[1];
rz(2.798583) q[1];
x q[2];
rz(1.0697332) q[3];
sx q[3];
rz(-2.2483629) q[3];
sx q[3];
rz(2.5354405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.23218368) q[2];
sx q[2];
rz(-2.3164618) q[2];
sx q[2];
rz(2.4278736) q[2];
rz(2.7632726) q[3];
sx q[3];
rz(-2.648073) q[3];
sx q[3];
rz(2.2617214) q[3];
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
rz(-pi) q[1];
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
rz(0.65080416) q[1];
sx q[1];
rz(-1.6497859) q[1];
sx q[1];
rz(-0.12129687) q[1];
rz(-1.7667608) q[2];
sx q[2];
rz(-1.6008196) q[2];
sx q[2];
rz(0.65884789) q[2];
rz(-1.1013423) q[3];
sx q[3];
rz(-2.6228842) q[3];
sx q[3];
rz(1.7026671) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];