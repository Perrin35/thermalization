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
rz(-0.38464889) q[0];
sx q[0];
rz(-2.3031213) q[0];
sx q[0];
rz(-2.5323001) q[0];
rz(0.72120178) q[1];
sx q[1];
rz(5.3546049) q[1];
sx q[1];
rz(7.3794853) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23966889) q[0];
sx q[0];
rz(-2.2623908) q[0];
sx q[0];
rz(2.9953629) q[0];
rz(-2.2112261) q[2];
sx q[2];
rz(-1.9631248) q[2];
sx q[2];
rz(1.0898255) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9723477) q[1];
sx q[1];
rz(-1.8512639) q[1];
sx q[1];
rz(-0.53963668) q[1];
rz(-pi) q[2];
rz(2.5240229) q[3];
sx q[3];
rz(-2.9520116) q[3];
sx q[3];
rz(-1.1763587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2234552) q[2];
sx q[2];
rz(-1.8310941) q[2];
sx q[2];
rz(1.1892595) q[2];
rz(0.39092815) q[3];
sx q[3];
rz(-1.1632185) q[3];
sx q[3];
rz(2.0093567) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44553462) q[0];
sx q[0];
rz(-1.3548509) q[0];
sx q[0];
rz(-1.498244) q[0];
rz(2.4636726) q[1];
sx q[1];
rz(-1.3005715) q[1];
sx q[1];
rz(-0.71151412) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2145793) q[0];
sx q[0];
rz(-1.2138146) q[0];
sx q[0];
rz(-2.0680319) q[0];
x q[1];
rz(-1.8078126) q[2];
sx q[2];
rz(-2.0227602) q[2];
sx q[2];
rz(2.1353561) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.38120684) q[1];
sx q[1];
rz(-2.7863656) q[1];
sx q[1];
rz(2.9095501) q[1];
rz(-0.031458843) q[3];
sx q[3];
rz(-0.94740564) q[3];
sx q[3];
rz(1.7707847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.70637643) q[2];
sx q[2];
rz(-1.7855568) q[2];
sx q[2];
rz(-1.0630652) q[2];
rz(1.6478018) q[3];
sx q[3];
rz(-2.6726275) q[3];
sx q[3];
rz(-2.4013605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4110334) q[0];
sx q[0];
rz(-2.7506802) q[0];
sx q[0];
rz(0.60182369) q[0];
rz(0.14566323) q[1];
sx q[1];
rz(-1.0258976) q[1];
sx q[1];
rz(-0.62785968) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3224761) q[0];
sx q[0];
rz(-1.6401344) q[0];
sx q[0];
rz(-0.11868166) q[0];
rz(1.026903) q[2];
sx q[2];
rz(-0.62335194) q[2];
sx q[2];
rz(-2.198213) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9722189) q[1];
sx q[1];
rz(-1.7986606) q[1];
sx q[1];
rz(-1.2078753) q[1];
x q[2];
rz(0.46247356) q[3];
sx q[3];
rz(-1.00849) q[3];
sx q[3];
rz(0.019788086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.38498983) q[2];
sx q[2];
rz(-1.0641229) q[2];
sx q[2];
rz(2.9260054) q[2];
rz(2.0032517) q[3];
sx q[3];
rz(-1.8214858) q[3];
sx q[3];
rz(-0.57083541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0067921) q[0];
sx q[0];
rz(-1.4147867) q[0];
sx q[0];
rz(-0.11882812) q[0];
rz(1.4745332) q[1];
sx q[1];
rz(-0.88913616) q[1];
sx q[1];
rz(-2.3337591) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0197163) q[0];
sx q[0];
rz(-2.3000679) q[0];
sx q[0];
rz(-0.76899935) q[0];
rz(-pi) q[1];
rz(3.0962871) q[2];
sx q[2];
rz(-1.2477562) q[2];
sx q[2];
rz(1.3458061) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4588759) q[1];
sx q[1];
rz(-1.7380875) q[1];
sx q[1];
rz(2.1840827) q[1];
x q[2];
rz(-1.586257) q[3];
sx q[3];
rz(-0.8925775) q[3];
sx q[3];
rz(-0.34567269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.64041758) q[2];
sx q[2];
rz(-1.7265604) q[2];
sx q[2];
rz(2.625722) q[2];
rz(0.094203146) q[3];
sx q[3];
rz(-2.3661864) q[3];
sx q[3];
rz(1.5532106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.77189389) q[0];
sx q[0];
rz(-2.0409245) q[0];
sx q[0];
rz(1.1871673) q[0];
rz(2.5502603) q[1];
sx q[1];
rz(-1.2966803) q[1];
sx q[1];
rz(-2.3147413) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2255838) q[0];
sx q[0];
rz(-1.1904612) q[0];
sx q[0];
rz(-1.7011786) q[0];
rz(-pi) q[1];
rz(3.060512) q[2];
sx q[2];
rz(-1.0141918) q[2];
sx q[2];
rz(-0.48559819) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3278153) q[1];
sx q[1];
rz(-2.2499488) q[1];
sx q[1];
rz(2.9742254) q[1];
rz(-pi) q[2];
rz(2.2777745) q[3];
sx q[3];
rz(-0.70138273) q[3];
sx q[3];
rz(0.4566628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.4225215) q[2];
sx q[2];
rz(-0.25455385) q[2];
sx q[2];
rz(1.0151601) q[2];
rz(2.2319131) q[3];
sx q[3];
rz(-1.2404975) q[3];
sx q[3];
rz(-0.66143405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61843094) q[0];
sx q[0];
rz(-1.2105415) q[0];
sx q[0];
rz(-2.8316408) q[0];
rz(0.018772086) q[1];
sx q[1];
rz(-1.680178) q[1];
sx q[1];
rz(-3.054256) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76742889) q[0];
sx q[0];
rz(-1.5383676) q[0];
sx q[0];
rz(-0.57291605) q[0];
x q[1];
rz(2.3051065) q[2];
sx q[2];
rz(-0.84794551) q[2];
sx q[2];
rz(1.370887) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0897905) q[1];
sx q[1];
rz(-0.6615122) q[1];
sx q[1];
rz(0.89349414) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0424006) q[3];
sx q[3];
rz(-1.434552) q[3];
sx q[3];
rz(0.50627499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.35149082) q[2];
sx q[2];
rz(-0.47372207) q[2];
sx q[2];
rz(2.6215485) q[2];
rz(-0.13488787) q[3];
sx q[3];
rz(-1.8205732) q[3];
sx q[3];
rz(-1.7322056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0670369) q[0];
sx q[0];
rz(-1.073607) q[0];
sx q[0];
rz(-2.7556162) q[0];
rz(1.0999673) q[1];
sx q[1];
rz(-2.4620582) q[1];
sx q[1];
rz(-2.3540672) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3271017) q[0];
sx q[0];
rz(-1.9270032) q[0];
sx q[0];
rz(0.30762958) q[0];
rz(-1.196127) q[2];
sx q[2];
rz(-0.9866937) q[2];
sx q[2];
rz(-1.0874334) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0213607) q[1];
sx q[1];
rz(-1.225705) q[1];
sx q[1];
rz(-2.3335275) q[1];
x q[2];
rz(0.16997108) q[3];
sx q[3];
rz(-1.8026226) q[3];
sx q[3];
rz(1.1500203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1197352) q[2];
sx q[2];
rz(-1.3081552) q[2];
sx q[2];
rz(1.5137399) q[2];
rz(-1.2210023) q[3];
sx q[3];
rz(-1.9767714) q[3];
sx q[3];
rz(2.6695719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8319703) q[0];
sx q[0];
rz(-1.4291052) q[0];
sx q[0];
rz(0.23304644) q[0];
rz(-1.4876935) q[1];
sx q[1];
rz(-1.033604) q[1];
sx q[1];
rz(1.0071365) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45389913) q[0];
sx q[0];
rz(-2.0845883) q[0];
sx q[0];
rz(-1.9625241) q[0];
rz(-2.51153) q[2];
sx q[2];
rz(-1.328397) q[2];
sx q[2];
rz(0.60464786) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5598232) q[1];
sx q[1];
rz(-1.1743465) q[1];
sx q[1];
rz(-0.15869402) q[1];
rz(-1.5910025) q[3];
sx q[3];
rz(-1.3267702) q[3];
sx q[3];
rz(-3.1283825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.86203376) q[2];
sx q[2];
rz(-2.0290012) q[2];
sx q[2];
rz(-0.88279185) q[2];
rz(0.27859303) q[3];
sx q[3];
rz(-1.168058) q[3];
sx q[3];
rz(2.6826503) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19926628) q[0];
sx q[0];
rz(-2.6658604) q[0];
sx q[0];
rz(1.5698154) q[0];
rz(0.78872952) q[1];
sx q[1];
rz(-2.1756344) q[1];
sx q[1];
rz(-0.68005651) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2325033) q[0];
sx q[0];
rz(-1.5111369) q[0];
sx q[0];
rz(0.76086126) q[0];
rz(-pi) q[1];
rz(1.6974738) q[2];
sx q[2];
rz(-1.4529723) q[2];
sx q[2];
rz(0.5165259) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9788821) q[1];
sx q[1];
rz(-1.0274402) q[1];
sx q[1];
rz(-1.1239978) q[1];
x q[2];
rz(0.96931501) q[3];
sx q[3];
rz(-0.82348862) q[3];
sx q[3];
rz(-1.4692509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1532229) q[2];
sx q[2];
rz(-0.92066568) q[2];
sx q[2];
rz(1.9752768) q[2];
rz(1.9370646) q[3];
sx q[3];
rz(-1.322999) q[3];
sx q[3];
rz(0.62674633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.523943) q[0];
sx q[0];
rz(-2.2594422) q[0];
sx q[0];
rz(2.7045265) q[0];
rz(0.79824671) q[1];
sx q[1];
rz(-0.50004807) q[1];
sx q[1];
rz(2.9214568) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94486124) q[0];
sx q[0];
rz(-1.5240977) q[0];
sx q[0];
rz(2.3406505) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2557477) q[2];
sx q[2];
rz(-1.3761259) q[2];
sx q[2];
rz(2.6597629) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70500532) q[1];
sx q[1];
rz(-1.844075) q[1];
sx q[1];
rz(2.336034) q[1];
rz(-pi) q[2];
rz(3.0978904) q[3];
sx q[3];
rz(-0.80207295) q[3];
sx q[3];
rz(1.5278096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.87307125) q[2];
sx q[2];
rz(-1.2348509) q[2];
sx q[2];
rz(0.29042563) q[2];
rz(-2.5051266) q[3];
sx q[3];
rz(-1.7593743) q[3];
sx q[3];
rz(1.2344454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81454043) q[0];
sx q[0];
rz(-2.4335813) q[0];
sx q[0];
rz(-2.0306564) q[0];
rz(-3.0771599) q[1];
sx q[1];
rz(-1.4705407) q[1];
sx q[1];
rz(2.4992117) q[1];
rz(-3.1101075) q[2];
sx q[2];
rz(-0.97618547) q[2];
sx q[2];
rz(0.57913274) q[2];
rz(-2.6024184) q[3];
sx q[3];
rz(-1.3414345) q[3];
sx q[3];
rz(0.10673005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
