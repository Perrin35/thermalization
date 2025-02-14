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
rz(0.60929259) q[0];
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
rz(-0.87920183) q[0];
sx q[0];
rz(-2.9953629) q[0];
rz(-pi) q[1];
rz(0.9303665) q[2];
sx q[2];
rz(-1.1784679) q[2];
sx q[2];
rz(2.0517672) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3071202) q[1];
sx q[1];
rz(-0.60170071) q[1];
sx q[1];
rz(2.6306319) q[1];
rz(-1.4601379) q[3];
sx q[3];
rz(-1.7250463) q[3];
sx q[3];
rz(-1.8024886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.91813749) q[2];
sx q[2];
rz(-1.8310941) q[2];
sx q[2];
rz(1.9523331) q[2];
rz(-2.7506645) q[3];
sx q[3];
rz(-1.1632185) q[3];
sx q[3];
rz(2.0093567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44553462) q[0];
sx q[0];
rz(-1.3548509) q[0];
sx q[0];
rz(-1.6433486) q[0];
rz(0.6779201) q[1];
sx q[1];
rz(-1.8410212) q[1];
sx q[1];
rz(2.4300785) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2133117) q[0];
sx q[0];
rz(-0.60325256) q[0];
sx q[0];
rz(2.2343982) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46320148) q[2];
sx q[2];
rz(-1.7836264) q[2];
sx q[2];
rz(-2.6821399) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.62809901) q[1];
sx q[1];
rz(-1.9160992) q[1];
sx q[1];
rz(-1.4856935) q[1];
rz(-pi) q[2];
rz(2.1944216) q[3];
sx q[3];
rz(-1.5452562) q[3];
sx q[3];
rz(-2.9232349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.70637643) q[2];
sx q[2];
rz(-1.3560359) q[2];
sx q[2];
rz(1.0630652) q[2];
rz(1.6478018) q[3];
sx q[3];
rz(-0.46896514) q[3];
sx q[3];
rz(-0.74023214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73055926) q[0];
sx q[0];
rz(-2.7506802) q[0];
sx q[0];
rz(-0.60182369) q[0];
rz(-0.14566323) q[1];
sx q[1];
rz(-2.115695) q[1];
sx q[1];
rz(2.513733) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8191166) q[0];
sx q[0];
rz(-1.5014582) q[0];
sx q[0];
rz(3.022911) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0192515) q[2];
sx q[2];
rz(-1.8776692) q[2];
sx q[2];
rz(-2.9705974) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9722189) q[1];
sx q[1];
rz(-1.7986606) q[1];
sx q[1];
rz(1.9337173) q[1];
x q[2];
rz(-0.95471621) q[3];
sx q[3];
rz(-0.71185116) q[3];
sx q[3];
rz(-2.4095132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.38498983) q[2];
sx q[2];
rz(-2.0774697) q[2];
sx q[2];
rz(0.21558726) q[2];
rz(1.138341) q[3];
sx q[3];
rz(-1.8214858) q[3];
sx q[3];
rz(0.57083541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1348006) q[0];
sx q[0];
rz(-1.726806) q[0];
sx q[0];
rz(0.11882812) q[0];
rz(1.6670594) q[1];
sx q[1];
rz(-2.2524565) q[1];
sx q[1];
rz(-2.3337591) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2654003) q[0];
sx q[0];
rz(-2.1159439) q[0];
sx q[0];
rz(-2.4643174) q[0];
rz(3.0962871) q[2];
sx q[2];
rz(-1.8938365) q[2];
sx q[2];
rz(1.7957866) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4588759) q[1];
sx q[1];
rz(-1.4035051) q[1];
sx q[1];
rz(2.1840827) q[1];
rz(-1.5553357) q[3];
sx q[3];
rz(-2.2490152) q[3];
sx q[3];
rz(2.79592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64041758) q[2];
sx q[2];
rz(-1.4150323) q[2];
sx q[2];
rz(0.51587063) q[2];
rz(0.094203146) q[3];
sx q[3];
rz(-0.7754063) q[3];
sx q[3];
rz(-1.5532106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.3696988) q[0];
sx q[0];
rz(-2.0409245) q[0];
sx q[0];
rz(1.1871673) q[0];
rz(2.5502603) q[1];
sx q[1];
rz(-1.8449123) q[1];
sx q[1];
rz(-0.82685131) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39385228) q[0];
sx q[0];
rz(-1.6918139) q[0];
sx q[0];
rz(0.38328247) q[0];
rz(-1.7002326) q[2];
sx q[2];
rz(-2.5797306) q[2];
sx q[2];
rz(2.8086086) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8137774) q[1];
sx q[1];
rz(-2.2499488) q[1];
sx q[1];
rz(-2.9742254) q[1];
x q[2];
rz(2.639797) q[3];
sx q[3];
rz(-1.0579946) q[3];
sx q[3];
rz(-0.38450188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7190711) q[2];
sx q[2];
rz(-0.25455385) q[2];
sx q[2];
rz(2.1264326) q[2];
rz(2.2319131) q[3];
sx q[3];
rz(-1.9010952) q[3];
sx q[3];
rz(0.66143405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61843094) q[0];
sx q[0];
rz(-1.9310512) q[0];
sx q[0];
rz(-0.30995187) q[0];
rz(3.1228206) q[1];
sx q[1];
rz(-1.4614146) q[1];
sx q[1];
rz(-3.054256) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3884443) q[0];
sx q[0];
rz(-0.57373057) q[0];
sx q[0];
rz(0.059772003) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87126731) q[2];
sx q[2];
rz(-2.0972898) q[2];
sx q[2];
rz(0.73824182) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8467056) q[1];
sx q[1];
rz(-1.0716096) q[1];
sx q[1];
rz(2.687672) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8363073) q[3];
sx q[3];
rz(-0.54406057) q[3];
sx q[3];
rz(1.8484704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7901018) q[2];
sx q[2];
rz(-2.6678706) q[2];
sx q[2];
rz(2.6215485) q[2];
rz(3.0067048) q[3];
sx q[3];
rz(-1.3210195) q[3];
sx q[3];
rz(-1.409387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0745558) q[0];
sx q[0];
rz(-1.073607) q[0];
sx q[0];
rz(-0.38597646) q[0];
rz(-1.0999673) q[1];
sx q[1];
rz(-2.4620582) q[1];
sx q[1];
rz(-0.78752548) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7875536) q[0];
sx q[0];
rz(-1.2830495) q[0];
sx q[0];
rz(-1.1985874) q[0];
rz(2.5239713) q[2];
sx q[2];
rz(-1.8810399) q[2];
sx q[2];
rz(-0.26981416) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0037076) q[1];
sx q[1];
rz(-2.2785997) q[1];
sx q[1];
rz(2.6801609) q[1];
x q[2];
rz(2.9716216) q[3];
sx q[3];
rz(-1.8026226) q[3];
sx q[3];
rz(1.9915723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0218574) q[2];
sx q[2];
rz(-1.3081552) q[2];
sx q[2];
rz(-1.6278527) q[2];
rz(1.2210023) q[3];
sx q[3];
rz(-1.9767714) q[3];
sx q[3];
rz(0.47202078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8319703) q[0];
sx q[0];
rz(-1.7124875) q[0];
sx q[0];
rz(-2.9085462) q[0];
rz(1.4876935) q[1];
sx q[1];
rz(-1.033604) q[1];
sx q[1];
rz(2.1344562) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45389913) q[0];
sx q[0];
rz(-1.0570044) q[0];
sx q[0];
rz(1.1790685) q[0];
x q[1];
rz(0.63006261) q[2];
sx q[2];
rz(-1.8131957) q[2];
sx q[2];
rz(2.5369448) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1888652) q[1];
sx q[1];
rz(-2.7161274) q[1];
sx q[1];
rz(-1.2098321) q[1];
rz(0.24407401) q[3];
sx q[3];
rz(-1.5511889) q[3];
sx q[3];
rz(-1.5527035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2795589) q[2];
sx q[2];
rz(-1.1125914) q[2];
sx q[2];
rz(-0.88279185) q[2];
rz(0.27859303) q[3];
sx q[3];
rz(-1.9735347) q[3];
sx q[3];
rz(-2.6826503) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9423264) q[0];
sx q[0];
rz(-2.6658604) q[0];
sx q[0];
rz(1.5717773) q[0];
rz(0.78872952) q[1];
sx q[1];
rz(-2.1756344) q[1];
sx q[1];
rz(-0.68005651) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4231789) q[0];
sx q[0];
rz(-2.3299651) q[0];
sx q[0];
rz(1.4885097) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4441188) q[2];
sx q[2];
rz(-1.4529723) q[2];
sx q[2];
rz(2.6250668) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.16526651) q[1];
sx q[1];
rz(-1.191947) q[1];
sx q[1];
rz(-0.59013175) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84362883) q[3];
sx q[3];
rz(-1.9988201) q[3];
sx q[3];
rz(-2.6035863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1532229) q[2];
sx q[2];
rz(-0.92066568) q[2];
sx q[2];
rz(-1.9752768) q[2];
rz(1.9370646) q[3];
sx q[3];
rz(-1.8185936) q[3];
sx q[3];
rz(2.5148463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.523943) q[0];
sx q[0];
rz(-2.2594422) q[0];
sx q[0];
rz(2.7045265) q[0];
rz(2.3433459) q[1];
sx q[1];
rz(-0.50004807) q[1];
sx q[1];
rz(-2.9214568) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58071346) q[0];
sx q[0];
rz(-0.80199837) q[0];
sx q[0];
rz(0.064994354) q[0];
rz(-pi) q[1];
rz(1.0041084) q[2];
sx q[2];
rz(-2.7729642) q[2];
sx q[2];
rz(1.6246375) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0221631) q[1];
sx q[1];
rz(-2.3009752) q[1];
sx q[1];
rz(2.7709211) q[1];
rz(-pi) q[2];
rz(1.5256568) q[3];
sx q[3];
rz(-2.3718826) q[3];
sx q[3];
rz(-1.5509645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2685214) q[2];
sx q[2];
rz(-1.2348509) q[2];
sx q[2];
rz(0.29042563) q[2];
rz(0.63646603) q[3];
sx q[3];
rz(-1.3822184) q[3];
sx q[3];
rz(-1.2344454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3270522) q[0];
sx q[0];
rz(-2.4335813) q[0];
sx q[0];
rz(-2.0306564) q[0];
rz(-3.0771599) q[1];
sx q[1];
rz(-1.4705407) q[1];
sx q[1];
rz(2.4992117) q[1];
rz(0.031485166) q[2];
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
