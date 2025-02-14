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
rz(-0.91312042) q[0];
sx q[0];
rz(-1.0364391) q[0];
sx q[0];
rz(1.6057462) q[0];
rz(-2.4802471) q[1];
sx q[1];
rz(-1.9547434) q[1];
sx q[1];
rz(-0.43651906) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8494328) q[0];
sx q[0];
rz(-1.8884553) q[0];
sx q[0];
rz(-3.0374797) q[0];
x q[1];
rz(0.5646602) q[2];
sx q[2];
rz(-1.1810978) q[2];
sx q[2];
rz(2.1291358) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.7750262) q[1];
sx q[1];
rz(-2.4373551) q[1];
sx q[1];
rz(-1.6338867) q[1];
rz(-pi) q[2];
rz(-0.17702509) q[3];
sx q[3];
rz(-0.89460556) q[3];
sx q[3];
rz(0.9392304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6785757) q[2];
sx q[2];
rz(-2.7028694) q[2];
sx q[2];
rz(-0.92301816) q[2];
rz(3.0758514) q[3];
sx q[3];
rz(-1.8140503) q[3];
sx q[3];
rz(1.8008697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0098669212) q[0];
sx q[0];
rz(-2.0781524) q[0];
sx q[0];
rz(-3.104082) q[0];
rz(0.58049479) q[1];
sx q[1];
rz(-0.66941222) q[1];
sx q[1];
rz(2.4494749) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1624533) q[0];
sx q[0];
rz(-1.7458436) q[0];
sx q[0];
rz(-1.2412423) q[0];
rz(-pi) q[1];
rz(-1.0672928) q[2];
sx q[2];
rz(-1.8177336) q[2];
sx q[2];
rz(1.4114789) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2381535) q[1];
sx q[1];
rz(-2.3336556) q[1];
sx q[1];
rz(-0.72407667) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7440657) q[3];
sx q[3];
rz(-0.82635802) q[3];
sx q[3];
rz(0.92191523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6633501) q[2];
sx q[2];
rz(-1.9448091) q[2];
sx q[2];
rz(1.7459858) q[2];
rz(2.9546402) q[3];
sx q[3];
rz(-1.170265) q[3];
sx q[3];
rz(-0.54005867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4726987) q[0];
sx q[0];
rz(-2.5085594) q[0];
sx q[0];
rz(2.582666) q[0];
rz(2.3129297) q[1];
sx q[1];
rz(-1.5507973) q[1];
sx q[1];
rz(-0.25159803) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64564451) q[0];
sx q[0];
rz(-2.6584266) q[0];
sx q[0];
rz(-0.42074411) q[0];
rz(-pi) q[1];
rz(2.0992804) q[2];
sx q[2];
rz(-2.5229075) q[2];
sx q[2];
rz(1.2109735) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.66183749) q[1];
sx q[1];
rz(-2.2244456) q[1];
sx q[1];
rz(-0.33262555) q[1];
rz(0.76520701) q[3];
sx q[3];
rz(-0.83929378) q[3];
sx q[3];
rz(0.30680233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.33637968) q[2];
sx q[2];
rz(-0.42625913) q[2];
sx q[2];
rz(1.7097998) q[2];
rz(-0.21670565) q[3];
sx q[3];
rz(-1.5313989) q[3];
sx q[3];
rz(1.7823035) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25877229) q[0];
sx q[0];
rz(-0.50676218) q[0];
sx q[0];
rz(-1.2910507) q[0];
rz(-2.3123815) q[1];
sx q[1];
rz(-2.0334838) q[1];
sx q[1];
rz(-1.0145899) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6651763) q[0];
sx q[0];
rz(-1.0422857) q[0];
sx q[0];
rz(1.2885874) q[0];
x q[1];
rz(0.56351557) q[2];
sx q[2];
rz(-1.9013828) q[2];
sx q[2];
rz(0.78621582) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17550771) q[1];
sx q[1];
rz(-2.1455742) q[1];
sx q[1];
rz(-0.0055343363) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7144373) q[3];
sx q[3];
rz(-0.82089409) q[3];
sx q[3];
rz(3.0526171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1318876) q[2];
sx q[2];
rz(-1.2010801) q[2];
sx q[2];
rz(0.76060549) q[2];
rz(2.6724114) q[3];
sx q[3];
rz(-1.4696591) q[3];
sx q[3];
rz(2.40707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.0932662) q[0];
sx q[0];
rz(-0.51437298) q[0];
sx q[0];
rz(-2.560428) q[0];
rz(-1.5515074) q[1];
sx q[1];
rz(-0.64684144) q[1];
sx q[1];
rz(2.4466628) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3669156) q[0];
sx q[0];
rz(-2.7645676) q[0];
sx q[0];
rz(0.89311262) q[0];
x q[1];
rz(-2.8492691) q[2];
sx q[2];
rz(-1.4806804) q[2];
sx q[2];
rz(2.9499049) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.1174005) q[1];
sx q[1];
rz(-1.0276737) q[1];
sx q[1];
rz(-1.622011) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7534689) q[3];
sx q[3];
rz(-1.338306) q[3];
sx q[3];
rz(0.35216613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.27802262) q[2];
sx q[2];
rz(-1.5333971) q[2];
sx q[2];
rz(-2.8387873) q[2];
rz(1.4263724) q[3];
sx q[3];
rz(-2.7808166) q[3];
sx q[3];
rz(-0.58437955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.67702174) q[0];
sx q[0];
rz(-0.40957054) q[0];
sx q[0];
rz(0.69497481) q[0];
rz(-2.8547844) q[1];
sx q[1];
rz(-0.60568714) q[1];
sx q[1];
rz(-1.274775) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0114637) q[0];
sx q[0];
rz(-1.090836) q[0];
sx q[0];
rz(0.35188727) q[0];
rz(-pi) q[1];
rz(2.5996501) q[2];
sx q[2];
rz(-1.9795609) q[2];
sx q[2];
rz(-2.8084076) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4732237) q[1];
sx q[1];
rz(-0.76938564) q[1];
sx q[1];
rz(1.046087) q[1];
rz(2.5596095) q[3];
sx q[3];
rz(-1.3701539) q[3];
sx q[3];
rz(2.0394005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0810762) q[2];
sx q[2];
rz(-0.21848564) q[2];
sx q[2];
rz(2.1742353) q[2];
rz(-0.71990144) q[3];
sx q[3];
rz(-2.1981809) q[3];
sx q[3];
rz(-1.2758183) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8571781) q[0];
sx q[0];
rz(-0.030817742) q[0];
sx q[0];
rz(1.4469294) q[0];
rz(-2.6743496) q[1];
sx q[1];
rz(-1.8485565) q[1];
sx q[1];
rz(1.1955059) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.104983) q[0];
sx q[0];
rz(-0.65147841) q[0];
sx q[0];
rz(-0.81654136) q[0];
rz(-pi) q[1];
rz(1.906997) q[2];
sx q[2];
rz(-2.3578491) q[2];
sx q[2];
rz(2.9072493) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6536234) q[1];
sx q[1];
rz(-1.5088827) q[1];
sx q[1];
rz(-2.0701522) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2310813) q[3];
sx q[3];
rz(-3.1034443) q[3];
sx q[3];
rz(1.3919786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.82470992) q[2];
sx q[2];
rz(-0.83493817) q[2];
sx q[2];
rz(2.6173124) q[2];
rz(2.8892062) q[3];
sx q[3];
rz(-3.0350244) q[3];
sx q[3];
rz(0.061080385) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97279945) q[0];
sx q[0];
rz(-1.1058829) q[0];
sx q[0];
rz(-1.9148781) q[0];
rz(0.75617689) q[1];
sx q[1];
rz(-0.94804472) q[1];
sx q[1];
rz(-0.39047584) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9928241) q[0];
sx q[0];
rz(-2.618194) q[0];
sx q[0];
rz(2.5958909) q[0];
rz(-1.5323213) q[2];
sx q[2];
rz(-2.0145825) q[2];
sx q[2];
rz(2.2930068) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9041216) q[1];
sx q[1];
rz(-1.7842222) q[1];
sx q[1];
rz(1.2280812) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49440439) q[3];
sx q[3];
rz(-0.98894152) q[3];
sx q[3];
rz(2.0928008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8172424) q[2];
sx q[2];
rz(-2.4411026) q[2];
sx q[2];
rz(-1.3520757) q[2];
rz(-2.9663441) q[3];
sx q[3];
rz(-1.8027571) q[3];
sx q[3];
rz(1.9372743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7290203) q[0];
sx q[0];
rz(-0.46171284) q[0];
sx q[0];
rz(-0.66226688) q[0];
rz(-0.63016713) q[1];
sx q[1];
rz(-1.8616734) q[1];
sx q[1];
rz(2.9060649) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92162392) q[0];
sx q[0];
rz(-1.3597836) q[0];
sx q[0];
rz(-1.3930265) q[0];
x q[1];
rz(-2.5723402) q[2];
sx q[2];
rz(-1.6496837) q[2];
sx q[2];
rz(-1.9663262) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3429361) q[1];
sx q[1];
rz(-1.6321811) q[1];
sx q[1];
rz(3.117901) q[1];
x q[2];
rz(3.0074545) q[3];
sx q[3];
rz(-2.2629177) q[3];
sx q[3];
rz(-1.9798756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0528637) q[2];
sx q[2];
rz(-0.95053089) q[2];
sx q[2];
rz(2.5843184) q[2];
rz(-2.3142464) q[3];
sx q[3];
rz(-1.6297623) q[3];
sx q[3];
rz(1.6415589) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2130704) q[0];
sx q[0];
rz(-2.3445573) q[0];
sx q[0];
rz(0.56761566) q[0];
rz(0.40191832) q[1];
sx q[1];
rz(-2.4979976) q[1];
sx q[1];
rz(1.7793122) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3338254) q[0];
sx q[0];
rz(-0.58192116) q[0];
sx q[0];
rz(-2.3533217) q[0];
x q[1];
rz(-1.9270114) q[2];
sx q[2];
rz(-2.491174) q[2];
sx q[2];
rz(-2.7678571) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4835254) q[1];
sx q[1];
rz(-2.1999252) q[1];
sx q[1];
rz(-0.097112819) q[1];
rz(-pi) q[2];
rz(0.74257039) q[3];
sx q[3];
rz(-0.97108632) q[3];
sx q[3];
rz(-2.7903583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.42084971) q[2];
sx q[2];
rz(-1.2030615) q[2];
sx q[2];
rz(0.26025772) q[2];
rz(0.69089729) q[3];
sx q[3];
rz(-0.8553718) q[3];
sx q[3];
rz(-1.6183841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33547587) q[0];
sx q[0];
rz(-1.5609043) q[0];
sx q[0];
rz(1.5326395) q[0];
rz(-2.7083022) q[1];
sx q[1];
rz(-1.3229803) q[1];
sx q[1];
rz(1.9569474) q[1];
rz(0.92930195) q[2];
sx q[2];
rz(-2.1718696) q[2];
sx q[2];
rz(1.0739506) q[2];
rz(-0.13570774) q[3];
sx q[3];
rz(-2.3466023) q[3];
sx q[3];
rz(0.96011163) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
