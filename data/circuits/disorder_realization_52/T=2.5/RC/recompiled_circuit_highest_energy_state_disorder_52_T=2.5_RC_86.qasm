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
rz(1.1531416) q[0];
sx q[0];
rz(-0.81557953) q[0];
sx q[0];
rz(2.3834035) q[0];
rz(-0.012501333) q[1];
sx q[1];
rz(4.9795436) q[1];
sx q[1];
rz(10.995168) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85938293) q[0];
sx q[0];
rz(-2.4042685) q[0];
sx q[0];
rz(-2.5928549) q[0];
rz(-pi) q[1];
rz(-1.7292132) q[2];
sx q[2];
rz(-2.7374501) q[2];
sx q[2];
rz(0.86106419) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4911035) q[1];
sx q[1];
rz(-1.2105618) q[1];
sx q[1];
rz(0.061042518) q[1];
rz(-pi) q[2];
rz(1.38405) q[3];
sx q[3];
rz(-2.8272326) q[3];
sx q[3];
rz(-1.7160742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5300753) q[2];
sx q[2];
rz(-3.1334183) q[2];
sx q[2];
rz(0.44069904) q[2];
rz(-0.079744451) q[3];
sx q[3];
rz(-0.00010448797) q[3];
sx q[3];
rz(1.124148) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1752862) q[0];
sx q[0];
rz(-3.0847302) q[0];
sx q[0];
rz(-2.9735907) q[0];
rz(-0.02027823) q[1];
sx q[1];
rz(-2.8333277) q[1];
sx q[1];
rz(1.6049989) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5791721) q[0];
sx q[0];
rz(-1.354166) q[0];
sx q[0];
rz(0.24026339) q[0];
rz(-pi) q[1];
rz(-1.3290908) q[2];
sx q[2];
rz(-3.1311718) q[2];
sx q[2];
rz(-1.7872321) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3520842) q[1];
sx q[1];
rz(-1.5696973) q[1];
sx q[1];
rz(1.5661245) q[1];
x q[2];
rz(1.5617989) q[3];
sx q[3];
rz(-1.6130884) q[3];
sx q[3];
rz(2.472714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7961879) q[2];
sx q[2];
rz(-2.2296843) q[2];
sx q[2];
rz(-1.4076642) q[2];
rz(2.0965072) q[3];
sx q[3];
rz(-3.0920691) q[3];
sx q[3];
rz(-2.869587) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3097836) q[0];
sx q[0];
rz(-2.164916) q[0];
sx q[0];
rz(-0.56104863) q[0];
rz(-2.8647515) q[1];
sx q[1];
rz(-3.1287153) q[1];
sx q[1];
rz(-1.8337839) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55371743) q[0];
sx q[0];
rz(-0.29969117) q[0];
sx q[0];
rz(1.4315579) q[0];
x q[1];
rz(-0.0073892825) q[2];
sx q[2];
rz(-1.5707914) q[2];
sx q[2];
rz(2.6342692) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6456994) q[1];
sx q[1];
rz(-0.99826854) q[1];
sx q[1];
rz(0.069620274) q[1];
rz(-pi) q[2];
rz(2.2482292) q[3];
sx q[3];
rz(-1.2033389) q[3];
sx q[3];
rz(1.8830255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7967367) q[2];
sx q[2];
rz(-0.00011809706) q[2];
sx q[2];
rz(2.5522088) q[2];
rz(1.899259) q[3];
sx q[3];
rz(-0.012367736) q[3];
sx q[3];
rz(1.7813659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0068483343) q[0];
sx q[0];
rz(-0.51181781) q[0];
sx q[0];
rz(-1.3486598) q[0];
rz(3.1349365) q[1];
sx q[1];
rz(-1.321188) q[1];
sx q[1];
rz(3.1080918) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5101227) q[0];
sx q[0];
rz(-0.92486787) q[0];
sx q[0];
rz(2.9246246) q[0];
rz(-pi) q[1];
rz(3.1042023) q[2];
sx q[2];
rz(-0.11971902) q[2];
sx q[2];
rz(0.38116383) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9814947) q[1];
sx q[1];
rz(-1.3019053) q[1];
sx q[1];
rz(0.010543028) q[1];
rz(-pi) q[2];
rz(1.4064404) q[3];
sx q[3];
rz(-1.4504315) q[3];
sx q[3];
rz(-3.0422378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1091619) q[2];
sx q[2];
rz(-3.1350632) q[2];
sx q[2];
rz(-2.8429441) q[2];
rz(-1.8715035) q[3];
sx q[3];
rz(-3.1254369) q[3];
sx q[3];
rz(-0.0531918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3700767) q[0];
sx q[0];
rz(-1.5144441) q[0];
sx q[0];
rz(2.694743) q[0];
rz(-0.18613786) q[1];
sx q[1];
rz(-0.061807241) q[1];
sx q[1];
rz(-1.7284547) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22767775) q[0];
sx q[0];
rz(-2.9083038) q[0];
sx q[0];
rz(-1.7167709) q[0];
rz(-pi) q[1];
rz(-1.9930219) q[2];
sx q[2];
rz(-1.7686426) q[2];
sx q[2];
rz(0.6126709) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6525938) q[1];
sx q[1];
rz(-1.6203383) q[1];
sx q[1];
rz(0.045157305) q[1];
rz(-pi) q[2];
rz(-0.90425332) q[3];
sx q[3];
rz(-0.35247856) q[3];
sx q[3];
rz(2.3942687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.83429217) q[2];
sx q[2];
rz(-1.5895546) q[2];
sx q[2];
rz(-2.6429122) q[2];
rz(-0.57791609) q[3];
sx q[3];
rz(-2.6579865) q[3];
sx q[3];
rz(2.5448866) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1805098) q[0];
sx q[0];
rz(-2.0208277) q[0];
sx q[0];
rz(2.7951796) q[0];
rz(-0.60180426) q[1];
sx q[1];
rz(-1.5806942) q[1];
sx q[1];
rz(0.7535038) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8418056) q[0];
sx q[0];
rz(-2.8827169) q[0];
sx q[0];
rz(-2.4095834) q[0];
x q[1];
rz(3.0664938) q[2];
sx q[2];
rz(-1.4595928) q[2];
sx q[2];
rz(-0.62125833) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1931455) q[1];
sx q[1];
rz(-0.79953927) q[1];
sx q[1];
rz(-2.1667797) q[1];
rz(-pi) q[2];
rz(3.0493204) q[3];
sx q[3];
rz(-1.7331976) q[3];
sx q[3];
rz(0.55264651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5715013) q[2];
sx q[2];
rz(-3.1381021) q[2];
sx q[2];
rz(-1.6174779) q[2];
rz(-3.0213455) q[3];
sx q[3];
rz(-0.0032987981) q[3];
sx q[3];
rz(-0.53774589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26871249) q[0];
sx q[0];
rz(-2.217642) q[0];
sx q[0];
rz(-0.22126108) q[0];
rz(-1.4606754) q[1];
sx q[1];
rz(-0.9333846) q[1];
sx q[1];
rz(0.077300765) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.59931) q[0];
sx q[0];
rz(-0.041754384) q[0];
sx q[0];
rz(1.6058654) q[0];
x q[1];
rz(1.5796698) q[2];
sx q[2];
rz(-1.5755782) q[2];
sx q[2];
rz(-1.3591131) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.12880023) q[1];
sx q[1];
rz(-0.18928738) q[1];
sx q[1];
rz(2.7526593) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12760136) q[3];
sx q[3];
rz(-1.7831047) q[3];
sx q[3];
rz(-1.9636167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7920502) q[2];
sx q[2];
rz(-0.011186102) q[2];
sx q[2];
rz(-2.1816317) q[2];
rz(-2.8121484) q[3];
sx q[3];
rz(-0.0080778413) q[3];
sx q[3];
rz(2.2913057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9462117) q[0];
sx q[0];
rz(-2.52849) q[0];
sx q[0];
rz(-0.10928133) q[0];
rz(2.7682313) q[1];
sx q[1];
rz(-0.80972087) q[1];
sx q[1];
rz(1.2304617) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62723535) q[0];
sx q[0];
rz(-2.1684596) q[0];
sx q[0];
rz(2.3742832) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3804761) q[2];
sx q[2];
rz(-1.7645823) q[2];
sx q[2];
rz(-1.5507669) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2956063) q[1];
sx q[1];
rz(-1.5593658) q[1];
sx q[1];
rz(1.5039526) q[1];
x q[2];
rz(-0.30970311) q[3];
sx q[3];
rz(-0.94325698) q[3];
sx q[3];
rz(2.5554339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5750778) q[2];
sx q[2];
rz(-1.2357624) q[2];
sx q[2];
rz(1.8129978) q[2];
rz(-1.7447507) q[3];
sx q[3];
rz(-0.003740398) q[3];
sx q[3];
rz(-1.0218792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-3.0777271) q[0];
sx q[0];
rz(-1.6985748) q[0];
sx q[0];
rz(0.57300895) q[0];
rz(2.833448) q[1];
sx q[1];
rz(-2.7317218) q[1];
sx q[1];
rz(-1.0073957) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1290179) q[0];
sx q[0];
rz(-1.5652302) q[0];
sx q[0];
rz(1.4642503) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0976866) q[2];
sx q[2];
rz(-1.708235) q[2];
sx q[2];
rz(-0.034417987) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3904265) q[1];
sx q[1];
rz(-1.6114283) q[1];
sx q[1];
rz(-1.4875814) q[1];
rz(-pi) q[2];
rz(1.151772) q[3];
sx q[3];
rz(-3.1085827) q[3];
sx q[3];
rz(-2.9275683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.31743) q[2];
sx q[2];
rz(-2.514826) q[2];
sx q[2];
rz(-2.7516348) q[2];
rz(-0.069084875) q[3];
sx q[3];
rz(-3.1324813) q[3];
sx q[3];
rz(-0.33153427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1214509) q[0];
sx q[0];
rz(-2.3905601) q[0];
sx q[0];
rz(2.6556515) q[0];
rz(-2.2700229) q[1];
sx q[1];
rz(-1.3078682) q[1];
sx q[1];
rz(-1.4942687) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6914929) q[0];
sx q[0];
rz(-1.0376103) q[0];
sx q[0];
rz(1.2144809) q[0];
x q[1];
rz(3.0790555) q[2];
sx q[2];
rz(-0.61574575) q[2];
sx q[2];
rz(-0.0074904895) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.941266) q[1];
sx q[1];
rz(-1.2709874) q[1];
sx q[1];
rz(-2.7977562) q[1];
rz(3.028454) q[3];
sx q[3];
rz(-1.5290878) q[3];
sx q[3];
rz(2.0737518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.5747052) q[2];
sx q[2];
rz(-0.042782728) q[2];
sx q[2];
rz(-3.1097143) q[2];
rz(-2.3632862) q[3];
sx q[3];
rz(-3.1347771) q[3];
sx q[3];
rz(-2.8460898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7185709) q[0];
sx q[0];
rz(-1.5324677) q[0];
sx q[0];
rz(1.8146429) q[0];
rz(0.12693916) q[1];
sx q[1];
rz(-0.23902421) q[1];
sx q[1];
rz(0.21993266) q[1];
rz(1.610582) q[2];
sx q[2];
rz(-3.0018158) q[2];
sx q[2];
rz(-2.9374585) q[2];
rz(0.056679139) q[3];
sx q[3];
rz(-2.2446742) q[3];
sx q[3];
rz(-2.6745924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
