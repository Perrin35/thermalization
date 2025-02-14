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
rz(0.83847133) q[0];
sx q[0];
rz(5.6738927) q[0];
rz(0.72120178) q[1];
sx q[1];
rz(-0.92858044) q[1];
sx q[1];
rz(-2.0452926) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1288777) q[0];
sx q[0];
rz(-0.70438526) q[0];
sx q[0];
rz(-1.3966068) q[0];
rz(-pi) q[1];
rz(0.9651411) q[2];
sx q[2];
rz(-0.7363626) q[2];
sx q[2];
rz(0.006927329) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23725736) q[1];
sx q[1];
rz(-1.0544027) q[1];
sx q[1];
rz(-1.2468491) q[1];
x q[2];
rz(-0.15518409) q[3];
sx q[3];
rz(-1.461457) q[3];
sx q[3];
rz(0.21462277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.91813749) q[2];
sx q[2];
rz(-1.8310941) q[2];
sx q[2];
rz(1.9523331) q[2];
rz(-0.39092815) q[3];
sx q[3];
rz(-1.1632185) q[3];
sx q[3];
rz(1.1322359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.696058) q[0];
sx q[0];
rz(-1.7867418) q[0];
sx q[0];
rz(1.498244) q[0];
rz(0.6779201) q[1];
sx q[1];
rz(-1.8410212) q[1];
sx q[1];
rz(-0.71151412) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9728015) q[0];
sx q[0];
rz(-2.0341221) q[0];
sx q[0];
rz(0.40131779) q[0];
rz(2.6783912) q[2];
sx q[2];
rz(-1.7836264) q[2];
sx q[2];
rz(2.6821399) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.62809901) q[1];
sx q[1];
rz(-1.9160992) q[1];
sx q[1];
rz(1.6558992) q[1];
rz(-pi) q[2];
rz(1.6145124) q[3];
sx q[3];
rz(-2.5175142) q[3];
sx q[3];
rz(1.3169552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4352162) q[2];
sx q[2];
rz(-1.3560359) q[2];
sx q[2];
rz(1.0630652) q[2];
rz(1.4937909) q[3];
sx q[3];
rz(-2.6726275) q[3];
sx q[3];
rz(-0.74023214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73055926) q[0];
sx q[0];
rz(-0.39091245) q[0];
sx q[0];
rz(2.539769) q[0];
rz(2.9959294) q[1];
sx q[1];
rz(-2.115695) q[1];
sx q[1];
rz(2.513733) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8635628) q[0];
sx q[0];
rz(-3.0042227) q[0];
sx q[0];
rz(0.53047116) q[0];
rz(-pi) q[1];
rz(2.1146896) q[2];
sx q[2];
rz(-0.62335194) q[2];
sx q[2];
rz(2.198213) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9380629) q[1];
sx q[1];
rz(-2.7157686) q[1];
sx q[1];
rz(0.99217207) q[1];
rz(2.1868764) q[3];
sx q[3];
rz(-2.4297415) q[3];
sx q[3];
rz(2.4095132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.38498983) q[2];
sx q[2];
rz(-2.0774697) q[2];
sx q[2];
rz(0.21558726) q[2];
rz(2.0032517) q[3];
sx q[3];
rz(-1.8214858) q[3];
sx q[3];
rz(2.5707572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0067921) q[0];
sx q[0];
rz(-1.726806) q[0];
sx q[0];
rz(-0.11882812) q[0];
rz(-1.4745332) q[1];
sx q[1];
rz(-2.2524565) q[1];
sx q[1];
rz(0.80783358) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2654003) q[0];
sx q[0];
rz(-2.1159439) q[0];
sx q[0];
rz(-0.67727526) q[0];
rz(-3.0962871) q[2];
sx q[2];
rz(-1.2477562) q[2];
sx q[2];
rz(-1.3458061) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4588759) q[1];
sx q[1];
rz(-1.4035051) q[1];
sx q[1];
rz(2.1840827) q[1];
rz(-pi) q[2];
rz(-0.67827722) q[3];
sx q[3];
rz(-1.5587574) q[3];
sx q[3];
rz(1.9067684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5011751) q[2];
sx q[2];
rz(-1.7265604) q[2];
sx q[2];
rz(0.51587063) q[2];
rz(-0.094203146) q[3];
sx q[3];
rz(-0.7754063) q[3];
sx q[3];
rz(1.5532106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3696988) q[0];
sx q[0];
rz(-1.1006681) q[0];
sx q[0];
rz(-1.1871673) q[0];
rz(0.59133235) q[1];
sx q[1];
rz(-1.8449123) q[1];
sx q[1];
rz(-2.3147413) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39385228) q[0];
sx q[0];
rz(-1.6918139) q[0];
sx q[0];
rz(-2.7583102) q[0];
x q[1];
rz(-3.060512) q[2];
sx q[2];
rz(-1.0141918) q[2];
sx q[2];
rz(0.48559819) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3278153) q[1];
sx q[1];
rz(-2.2499488) q[1];
sx q[1];
rz(2.9742254) q[1];
x q[2];
rz(-2.1416792) q[3];
sx q[3];
rz(-1.1383082) q[3];
sx q[3];
rz(-1.4492346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7190711) q[2];
sx q[2];
rz(-0.25455385) q[2];
sx q[2];
rz(2.1264326) q[2];
rz(-2.2319131) q[3];
sx q[3];
rz(-1.9010952) q[3];
sx q[3];
rz(-0.66143405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61843094) q[0];
sx q[0];
rz(-1.2105415) q[0];
sx q[0];
rz(-0.30995187) q[0];
rz(0.018772086) q[1];
sx q[1];
rz(-1.680178) q[1];
sx q[1];
rz(0.087336691) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7531484) q[0];
sx q[0];
rz(-0.57373057) q[0];
sx q[0];
rz(0.059772003) q[0];
rz(2.3051065) q[2];
sx q[2];
rz(-2.2936471) q[2];
sx q[2];
rz(-1.370887) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0951335) q[1];
sx q[1];
rz(-1.1756056) q[1];
sx q[1];
rz(1.0254775) q[1];
rz(-1.0424006) q[3];
sx q[3];
rz(-1.434552) q[3];
sx q[3];
rz(-0.50627499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.35149082) q[2];
sx q[2];
rz(-2.6678706) q[2];
sx q[2];
rz(2.6215485) q[2];
rz(0.13488787) q[3];
sx q[3];
rz(-1.8205732) q[3];
sx q[3];
rz(1.7322056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0670369) q[0];
sx q[0];
rz(-1.073607) q[0];
sx q[0];
rz(-0.38597646) q[0];
rz(-2.0416253) q[1];
sx q[1];
rz(-0.67953449) q[1];
sx q[1];
rz(-0.78752548) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7875536) q[0];
sx q[0];
rz(-1.2830495) q[0];
sx q[0];
rz(-1.1985874) q[0];
rz(-pi) q[1];
rz(2.5239713) q[2];
sx q[2];
rz(-1.8810399) q[2];
sx q[2];
rz(-0.26981416) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.13788504) q[1];
sx q[1];
rz(-0.862993) q[1];
sx q[1];
rz(-2.6801609) q[1];
x q[2];
rz(0.94905628) q[3];
sx q[3];
rz(-2.8550386) q[3];
sx q[3];
rz(-1.3499945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0218574) q[2];
sx q[2];
rz(-1.8334374) q[2];
sx q[2];
rz(-1.6278527) q[2];
rz(-1.2210023) q[3];
sx q[3];
rz(-1.1648213) q[3];
sx q[3];
rz(0.47202078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8319703) q[0];
sx q[0];
rz(-1.7124875) q[0];
sx q[0];
rz(-2.9085462) q[0];
rz(1.6538992) q[1];
sx q[1];
rz(-1.033604) q[1];
sx q[1];
rz(1.0071365) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45389913) q[0];
sx q[0];
rz(-2.0845883) q[0];
sx q[0];
rz(1.1790685) q[0];
rz(-2.7442544) q[2];
sx q[2];
rz(-0.66911829) q[2];
sx q[2];
rz(2.4934409) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9527275) q[1];
sx q[1];
rz(-2.7161274) q[1];
sx q[1];
rz(1.2098321) q[1];
x q[2];
rz(-1.5910025) q[3];
sx q[3];
rz(-1.3267702) q[3];
sx q[3];
rz(0.013210162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2795589) q[2];
sx q[2];
rz(-2.0290012) q[2];
sx q[2];
rz(-2.2588008) q[2];
rz(-2.8629996) q[3];
sx q[3];
rz(-1.168058) q[3];
sx q[3];
rz(2.6826503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19926628) q[0];
sx q[0];
rz(-2.6658604) q[0];
sx q[0];
rz(1.5717773) q[0];
rz(2.3528631) q[1];
sx q[1];
rz(-2.1756344) q[1];
sx q[1];
rz(-2.4615361) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4231789) q[0];
sx q[0];
rz(-0.81162757) q[0];
sx q[0];
rz(1.4885097) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4441188) q[2];
sx q[2];
rz(-1.6886204) q[2];
sx q[2];
rz(2.6250668) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1627106) q[1];
sx q[1];
rz(-1.0274402) q[1];
sx q[1];
rz(-2.0175949) q[1];
rz(-pi) q[2];
rz(0.96931501) q[3];
sx q[3];
rz(-0.82348862) q[3];
sx q[3];
rz(1.6723417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1532229) q[2];
sx q[2];
rz(-0.92066568) q[2];
sx q[2];
rz(-1.1663158) q[2];
rz(-1.204528) q[3];
sx q[3];
rz(-1.322999) q[3];
sx q[3];
rz(-2.5148463) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.523943) q[0];
sx q[0];
rz(-0.88215041) q[0];
sx q[0];
rz(0.43706617) q[0];
rz(2.3433459) q[1];
sx q[1];
rz(-0.50004807) q[1];
sx q[1];
rz(0.22013586) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6740538) q[0];
sx q[0];
rz(-2.3706145) q[0];
sx q[0];
rz(-1.6378372) q[0];
rz(2.1374843) q[2];
sx q[2];
rz(-0.36862843) q[2];
sx q[2];
rz(1.6246375) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4365873) q[1];
sx q[1];
rz(-1.2975177) q[1];
sx q[1];
rz(0.80555861) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0978904) q[3];
sx q[3];
rz(-0.80207295) q[3];
sx q[3];
rz(1.5278096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.87307125) q[2];
sx q[2];
rz(-1.2348509) q[2];
sx q[2];
rz(2.851167) q[2];
rz(-2.5051266) q[3];
sx q[3];
rz(-1.3822184) q[3];
sx q[3];
rz(1.9071473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3270522) q[0];
sx q[0];
rz(-0.70801133) q[0];
sx q[0];
rz(1.1109362) q[0];
rz(0.064432714) q[1];
sx q[1];
rz(-1.4705407) q[1];
sx q[1];
rz(2.4992117) q[1];
rz(3.1101075) q[2];
sx q[2];
rz(-2.1654072) q[2];
sx q[2];
rz(-2.5624599) q[2];
rz(-0.42677943) q[3];
sx q[3];
rz(-2.5601294) q[3];
sx q[3];
rz(-1.1009205) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
