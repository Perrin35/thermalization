OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0107083) q[0];
sx q[0];
rz(-2.4680128) q[0];
sx q[0];
rz(-2.3186865) q[0];
rz(1.4505439) q[1];
sx q[1];
rz(-0.6757285) q[1];
sx q[1];
rz(-2.9339209) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9638478) q[0];
sx q[0];
rz(-1.9152993) q[0];
sx q[0];
rz(0.81015195) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0479716) q[2];
sx q[2];
rz(-2.516921) q[2];
sx q[2];
rz(-0.60909787) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.80245167) q[1];
sx q[1];
rz(-1.4467518) q[1];
sx q[1];
rz(-1.63271) q[1];
x q[2];
rz(-3.0700686) q[3];
sx q[3];
rz(-1.9404963) q[3];
sx q[3];
rz(0.92310753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1755918) q[2];
sx q[2];
rz(-1.5755743) q[2];
sx q[2];
rz(-1.9123745) q[2];
rz(1.843533) q[3];
sx q[3];
rz(-1.3763873) q[3];
sx q[3];
rz(-1.1840597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17756473) q[0];
sx q[0];
rz(-0.25968817) q[0];
sx q[0];
rz(-0.92726707) q[0];
rz(-2.941653) q[1];
sx q[1];
rz(-1.4727458) q[1];
sx q[1];
rz(0.98958579) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7829076) q[0];
sx q[0];
rz(-1.9316422) q[0];
sx q[0];
rz(-3.0263958) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.093002158) q[2];
sx q[2];
rz(-0.43605294) q[2];
sx q[2];
rz(-2.2047037) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.687508) q[1];
sx q[1];
rz(-1.6191779) q[1];
sx q[1];
rz(1.494074) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2842312) q[3];
sx q[3];
rz(-1.1097842) q[3];
sx q[3];
rz(1.1198514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2135311) q[2];
sx q[2];
rz(-1.7240883) q[2];
sx q[2];
rz(0.35378635) q[2];
rz(0.95156041) q[3];
sx q[3];
rz(-2.3944201) q[3];
sx q[3];
rz(-0.32683867) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3554409) q[0];
sx q[0];
rz(-0.94796258) q[0];
sx q[0];
rz(1.0107262) q[0];
rz(-0.058723681) q[1];
sx q[1];
rz(-1.6207691) q[1];
sx q[1];
rz(-0.40930632) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.302251) q[0];
sx q[0];
rz(-0.88345655) q[0];
sx q[0];
rz(0.99498596) q[0];
rz(-pi) q[1];
rz(1.204973) q[2];
sx q[2];
rz(-0.43667781) q[2];
sx q[2];
rz(0.15966378) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5262437) q[1];
sx q[1];
rz(-1.0305291) q[1];
sx q[1];
rz(1.8002285) q[1];
rz(-1.5915967) q[3];
sx q[3];
rz(-1.0619508) q[3];
sx q[3];
rz(-1.5326981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.056444082) q[2];
sx q[2];
rz(-1.2480382) q[2];
sx q[2];
rz(0.15963456) q[2];
rz(1.8917482) q[3];
sx q[3];
rz(-1.4188473) q[3];
sx q[3];
rz(1.0861402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98642629) q[0];
sx q[0];
rz(-2.7041589) q[0];
sx q[0];
rz(1.0945818) q[0];
rz(-1.1711586) q[1];
sx q[1];
rz(-2.0234225) q[1];
sx q[1];
rz(-1.0135244) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60166925) q[0];
sx q[0];
rz(-0.058383103) q[0];
sx q[0];
rz(0.71805449) q[0];
rz(-pi) q[1];
rz(-1.5640352) q[2];
sx q[2];
rz(-2.6471615) q[2];
sx q[2];
rz(1.0552366) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3588874) q[1];
sx q[1];
rz(-0.26391477) q[1];
sx q[1];
rz(0.15586075) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0059399) q[3];
sx q[3];
rz(-1.4732822) q[3];
sx q[3];
rz(1.2312191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0052428) q[2];
sx q[2];
rz(-1.3815657) q[2];
sx q[2];
rz(1.8278149) q[2];
rz(2.2037196) q[3];
sx q[3];
rz(-1.8504986) q[3];
sx q[3];
rz(3.042799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90792847) q[0];
sx q[0];
rz(-1.5282682) q[0];
sx q[0];
rz(2.0966356) q[0];
rz(-2.9772671) q[1];
sx q[1];
rz(-1.3427837) q[1];
sx q[1];
rz(2.9716861) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1327672) q[0];
sx q[0];
rz(-0.19951162) q[0];
sx q[0];
rz(1.2037174) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3218574) q[2];
sx q[2];
rz(-2.6728874) q[2];
sx q[2];
rz(-0.021989487) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2881238) q[1];
sx q[1];
rz(-0.92902029) q[1];
sx q[1];
rz(-0.30499129) q[1];
rz(1.4395797) q[3];
sx q[3];
rz(-0.84892143) q[3];
sx q[3];
rz(3.1027681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5728411) q[2];
sx q[2];
rz(-0.70256394) q[2];
sx q[2];
rz(-1.0315726) q[2];
rz(1.0860363) q[3];
sx q[3];
rz(-0.7998172) q[3];
sx q[3];
rz(1.4097376) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80424911) q[0];
sx q[0];
rz(-2.0954837) q[0];
sx q[0];
rz(-1.7403437) q[0];
rz(2.7119472) q[1];
sx q[1];
rz(-2.255217) q[1];
sx q[1];
rz(1.3053798) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3004634) q[0];
sx q[0];
rz(-0.80781898) q[0];
sx q[0];
rz(-0.43231583) q[0];
rz(-1.5804703) q[2];
sx q[2];
rz(-2.7026483) q[2];
sx q[2];
rz(-2.2679971) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3363477) q[1];
sx q[1];
rz(-1.3224396) q[1];
sx q[1];
rz(-1.6026366) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8135728) q[3];
sx q[3];
rz(-1.4221898) q[3];
sx q[3];
rz(0.44156238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1244916) q[2];
sx q[2];
rz(-1.2064826) q[2];
sx q[2];
rz(2.1567832) q[2];
rz(-1.56636) q[3];
sx q[3];
rz(-1.5896348) q[3];
sx q[3];
rz(0.28520939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.259909) q[0];
sx q[0];
rz(-0.30170983) q[0];
sx q[0];
rz(-1.3551706) q[0];
rz(-3.1175218) q[1];
sx q[1];
rz(-1.6203974) q[1];
sx q[1];
rz(0.37758652) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51272362) q[0];
sx q[0];
rz(-2.2170904) q[0];
sx q[0];
rz(1.3747526) q[0];
rz(0.88084282) q[2];
sx q[2];
rz(-0.27327785) q[2];
sx q[2];
rz(2.5246594) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7166669) q[1];
sx q[1];
rz(-1.9583988) q[1];
sx q[1];
rz(-2.7295154) q[1];
rz(-pi) q[2];
rz(-2.2144775) q[3];
sx q[3];
rz(-2.6171631) q[3];
sx q[3];
rz(-0.35628179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5923803) q[2];
sx q[2];
rz(-2.8040631) q[2];
sx q[2];
rz(0.40360061) q[2];
rz(0.18925439) q[3];
sx q[3];
rz(-1.0523825) q[3];
sx q[3];
rz(-2.6846867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6117578) q[0];
sx q[0];
rz(-2.5921322) q[0];
sx q[0];
rz(-2.8305565) q[0];
rz(-0.95651904) q[1];
sx q[1];
rz(-1.4776968) q[1];
sx q[1];
rz(0.32435736) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8984025) q[0];
sx q[0];
rz(-2.0721738) q[0];
sx q[0];
rz(-1.7883975) q[0];
rz(-pi) q[1];
rz(-3.0359984) q[2];
sx q[2];
rz(-2.6289231) q[2];
sx q[2];
rz(-1.6216506) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.94453401) q[1];
sx q[1];
rz(-2.1558216) q[1];
sx q[1];
rz(0.50307517) q[1];
rz(-pi) q[2];
rz(-2.0288386) q[3];
sx q[3];
rz(-0.73087091) q[3];
sx q[3];
rz(-0.79938408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.73584622) q[2];
sx q[2];
rz(-0.79068557) q[2];
sx q[2];
rz(2.9131367) q[2];
rz(2.6089846) q[3];
sx q[3];
rz(-2.3753128) q[3];
sx q[3];
rz(-2.2920091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5936977) q[0];
sx q[0];
rz(-0.51012796) q[0];
sx q[0];
rz(1.2702031) q[0];
rz(0.58865976) q[1];
sx q[1];
rz(-1.9344067) q[1];
sx q[1];
rz(2.7587845) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8981738) q[0];
sx q[0];
rz(-0.23056689) q[0];
sx q[0];
rz(-1.6797934) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50254681) q[2];
sx q[2];
rz(-2.655173) q[2];
sx q[2];
rz(2.9596299) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3119024) q[1];
sx q[1];
rz(-1.6964579) q[1];
sx q[1];
rz(2.4290393) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8715897) q[3];
sx q[3];
rz(-2.0473891) q[3];
sx q[3];
rz(-1.4780413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1560912) q[2];
sx q[2];
rz(-1.568855) q[2];
sx q[2];
rz(2.7413979) q[2];
rz(-0.24724809) q[3];
sx q[3];
rz(-1.4618382) q[3];
sx q[3];
rz(0.39321536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-1.8409214) q[0];
sx q[0];
rz(-1.8074169) q[0];
sx q[0];
rz(-2.1260496) q[0];
rz(-2.4726942) q[1];
sx q[1];
rz(-1.6872419) q[1];
sx q[1];
rz(1.5060172) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7114764) q[0];
sx q[0];
rz(-1.073494) q[0];
sx q[0];
rz(-0.74691746) q[0];
rz(-0.54803063) q[2];
sx q[2];
rz(-1.9491674) q[2];
sx q[2];
rz(0.038383287) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.64262154) q[1];
sx q[1];
rz(-1.5318512) q[1];
sx q[1];
rz(1.0907111) q[1];
rz(-pi) q[2];
rz(0.56699879) q[3];
sx q[3];
rz(-2.1666489) q[3];
sx q[3];
rz(-1.8351549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.16613913) q[2];
sx q[2];
rz(-2.0659122) q[2];
sx q[2];
rz(1.6157185) q[2];
rz(-2.8499991) q[3];
sx q[3];
rz(-2.6988131) q[3];
sx q[3];
rz(2.7899138) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6780554) q[0];
sx q[0];
rz(-2.1778477) q[0];
sx q[0];
rz(0.5974593) q[0];
rz(-0.28868227) q[1];
sx q[1];
rz(-2.3434227) q[1];
sx q[1];
rz(-3.1176288) q[1];
rz(-1.8815251) q[2];
sx q[2];
rz(-1.962553) q[2];
sx q[2];
rz(0.16061781) q[2];
rz(-1.3194094) q[3];
sx q[3];
rz(-1.1234) q[3];
sx q[3];
rz(3.1004536) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
