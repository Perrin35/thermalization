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
rz(-1.5795213) q[0];
sx q[0];
rz(-3.002394) q[0];
sx q[0];
rz(2.8237999) q[0];
rz(-2.3843482) q[1];
sx q[1];
rz(4.7321893) q[1];
sx q[1];
rz(7.7319747) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23148558) q[0];
sx q[0];
rz(-0.39662374) q[0];
sx q[0];
rz(-2.7825732) q[0];
rz(-1.4124539) q[2];
sx q[2];
rz(-1.2333561) q[2];
sx q[2];
rz(-0.98357633) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6283012) q[1];
sx q[1];
rz(-2.5963497) q[1];
sx q[1];
rz(-2.1907275) q[1];
rz(-pi) q[2];
rz(1.3643614) q[3];
sx q[3];
rz(-0.41461333) q[3];
sx q[3];
rz(1.9569091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.455787) q[2];
sx q[2];
rz(-2.0152342) q[2];
sx q[2];
rz(0.2618928) q[2];
rz(0.68771466) q[3];
sx q[3];
rz(-2.3759638) q[3];
sx q[3];
rz(-0.31622893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19551109) q[0];
sx q[0];
rz(-1.8571778) q[0];
sx q[0];
rz(2.2971357) q[0];
rz(0.63956815) q[1];
sx q[1];
rz(-0.25194672) q[1];
sx q[1];
rz(-1.8738939) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9828192) q[0];
sx q[0];
rz(-0.70342677) q[0];
sx q[0];
rz(1.7577459) q[0];
rz(-pi) q[1];
rz(0.19494178) q[2];
sx q[2];
rz(-0.6856853) q[2];
sx q[2];
rz(2.1491988) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4417394) q[1];
sx q[1];
rz(-1.5754712) q[1];
sx q[1];
rz(2.38928) q[1];
rz(0.58749871) q[3];
sx q[3];
rz(-2.2031257) q[3];
sx q[3];
rz(-1.895623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1203479) q[2];
sx q[2];
rz(-0.52560386) q[2];
sx q[2];
rz(0.28908238) q[2];
rz(-2.8287079) q[3];
sx q[3];
rz(-0.94066921) q[3];
sx q[3];
rz(0.59278929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4631571) q[0];
sx q[0];
rz(-2.1364697) q[0];
sx q[0];
rz(-2.4682755) q[0];
rz(2.8939269) q[1];
sx q[1];
rz(-2.1238056) q[1];
sx q[1];
rz(0.30390513) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9760151) q[0];
sx q[0];
rz(-1.0212693) q[0];
sx q[0];
rz(2.5499792) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2938305) q[2];
sx q[2];
rz(-2.1925732) q[2];
sx q[2];
rz(1.997681) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.9527377) q[1];
sx q[1];
rz(-1.6293289) q[1];
sx q[1];
rz(-1.7105647) q[1];
rz(-0.21411333) q[3];
sx q[3];
rz(-1.1592093) q[3];
sx q[3];
rz(2.0579091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5725382) q[2];
sx q[2];
rz(-1.7819541) q[2];
sx q[2];
rz(-3.1247826) q[2];
rz(2.861764) q[3];
sx q[3];
rz(-0.40585104) q[3];
sx q[3];
rz(2.5376937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0458321) q[0];
sx q[0];
rz(-2.026676) q[0];
sx q[0];
rz(1.8810077) q[0];
rz(-0.49651217) q[1];
sx q[1];
rz(-1.9910201) q[1];
sx q[1];
rz(1.6373434) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53816089) q[0];
sx q[0];
rz(-1.0350845) q[0];
sx q[0];
rz(-2.6441022) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8102747) q[2];
sx q[2];
rz(-2.1000855) q[2];
sx q[2];
rz(-0.052241922) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.49676412) q[1];
sx q[1];
rz(-1.033492) q[1];
sx q[1];
rz(0.19283466) q[1];
rz(-pi) q[2];
rz(-2.6332985) q[3];
sx q[3];
rz(-2.5807305) q[3];
sx q[3];
rz(1.1713374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0513976) q[2];
sx q[2];
rz(-0.53956705) q[2];
sx q[2];
rz(0.16057333) q[2];
rz(-1.4422902) q[3];
sx q[3];
rz(-1.620404) q[3];
sx q[3];
rz(-0.035339385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(0.2546688) q[0];
sx q[0];
rz(-1.2441607) q[0];
sx q[0];
rz(-0.96920335) q[0];
rz(1.8907549) q[1];
sx q[1];
rz(-2.4617317) q[1];
sx q[1];
rz(0.23076898) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7353775) q[0];
sx q[0];
rz(-3.1054524) q[0];
sx q[0];
rz(-1.7202431) q[0];
rz(-1.5838301) q[2];
sx q[2];
rz(-0.37838867) q[2];
sx q[2];
rz(-0.8410483) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.90143067) q[1];
sx q[1];
rz(-2.0367176) q[1];
sx q[1];
rz(0.45216143) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67471116) q[3];
sx q[3];
rz(-2.5541948) q[3];
sx q[3];
rz(0.34492465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1988924) q[2];
sx q[2];
rz(-0.57622272) q[2];
sx q[2];
rz(2.8885941) q[2];
rz(1.7871855) q[3];
sx q[3];
rz(-1.1043045) q[3];
sx q[3];
rz(-1.3033006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81786466) q[0];
sx q[0];
rz(-2.4619894) q[0];
sx q[0];
rz(3.0815109) q[0];
rz(0.59728638) q[1];
sx q[1];
rz(-2.2164454) q[1];
sx q[1];
rz(3.0937321) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40210783) q[0];
sx q[0];
rz(-0.093221752) q[0];
sx q[0];
rz(-1.4076821) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5563407) q[2];
sx q[2];
rz(-1.774228) q[2];
sx q[2];
rz(-2.360441) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7737433) q[1];
sx q[1];
rz(-1.7715104) q[1];
sx q[1];
rz(-1.9108921) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32190692) q[3];
sx q[3];
rz(-1.6471146) q[3];
sx q[3];
rz(0.25086153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4045599) q[2];
sx q[2];
rz(-0.76475778) q[2];
sx q[2];
rz(0.16727373) q[2];
rz(-0.078992756) q[3];
sx q[3];
rz(-1.9588214) q[3];
sx q[3];
rz(2.3110068) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1792184) q[0];
sx q[0];
rz(-0.11070624) q[0];
sx q[0];
rz(0.12744823) q[0];
rz(2.2178862) q[1];
sx q[1];
rz(-0.5674924) q[1];
sx q[1];
rz(-2.4097402) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6293242) q[0];
sx q[0];
rz(-1.4254117) q[0];
sx q[0];
rz(-0.54319197) q[0];
rz(-1.8867854) q[2];
sx q[2];
rz(-1.1605186) q[2];
sx q[2];
rz(2.7261594) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4036847) q[1];
sx q[1];
rz(-2.2353454) q[1];
sx q[1];
rz(0.43364201) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50639373) q[3];
sx q[3];
rz(-0.78575883) q[3];
sx q[3];
rz(2.9036811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7544516) q[2];
sx q[2];
rz(-1.0820505) q[2];
sx q[2];
rz(2.2608536) q[2];
rz(-2.8494917) q[3];
sx q[3];
rz(-0.62551522) q[3];
sx q[3];
rz(2.1706799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3683415) q[0];
sx q[0];
rz(-2.1357949) q[0];
sx q[0];
rz(0.7811501) q[0];
rz(-1.7225522) q[1];
sx q[1];
rz(-2.1731989) q[1];
sx q[1];
rz(0.17793812) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74465655) q[0];
sx q[0];
rz(-1.2632003) q[0];
sx q[0];
rz(1.164695) q[0];
rz(-pi) q[1];
rz(-0.85748144) q[2];
sx q[2];
rz(-2.4861587) q[2];
sx q[2];
rz(-2.0542184) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3314245) q[1];
sx q[1];
rz(-1.6292062) q[1];
sx q[1];
rz(2.057079) q[1];
x q[2];
rz(-0.26375651) q[3];
sx q[3];
rz(-1.7510771) q[3];
sx q[3];
rz(-0.83747415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.51131821) q[2];
sx q[2];
rz(-1.7061468) q[2];
sx q[2];
rz(1.1159631) q[2];
rz(2.337194) q[3];
sx q[3];
rz(-2.4837327) q[3];
sx q[3];
rz(-0.72949725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6804752) q[0];
sx q[0];
rz(-0.68515968) q[0];
sx q[0];
rz(-2.6035736) q[0];
rz(2.8021725) q[1];
sx q[1];
rz(-1.5433886) q[1];
sx q[1];
rz(3.0982049) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8755084) q[0];
sx q[0];
rz(-1.5055297) q[0];
sx q[0];
rz(1.6677744) q[0];
rz(-pi) q[1];
x q[1];
rz(1.814079) q[2];
sx q[2];
rz(-2.5777811) q[2];
sx q[2];
rz(2.3772637) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.63396066) q[1];
sx q[1];
rz(-1.6502737) q[1];
sx q[1];
rz(0.75130557) q[1];
rz(-2.6049117) q[3];
sx q[3];
rz(-2.4701475) q[3];
sx q[3];
rz(2.1083903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.033120774) q[2];
sx q[2];
rz(-0.35999173) q[2];
sx q[2];
rz(1.6610422) q[2];
rz(0.2964274) q[3];
sx q[3];
rz(-1.1400283) q[3];
sx q[3];
rz(-2.3886783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4794889) q[0];
sx q[0];
rz(-0.39283735) q[0];
sx q[0];
rz(1.1641077) q[0];
rz(-0.33319831) q[1];
sx q[1];
rz(-1.6644128) q[1];
sx q[1];
rz(2.3936757) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18834297) q[0];
sx q[0];
rz(-0.97916257) q[0];
sx q[0];
rz(1.9337898) q[0];
rz(-pi) q[1];
rz(-2.9299253) q[2];
sx q[2];
rz(-0.94918409) q[2];
sx q[2];
rz(1.399161) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.235215) q[1];
sx q[1];
rz(-3.0864419) q[1];
sx q[1];
rz(-0.13184594) q[1];
rz(-1.6259057) q[3];
sx q[3];
rz(-0.69384533) q[3];
sx q[3];
rz(1.2220647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.36772874) q[2];
sx q[2];
rz(-1.7325956) q[2];
sx q[2];
rz(-1.1665974) q[2];
rz(-1.599132) q[3];
sx q[3];
rz(-1.6049933) q[3];
sx q[3];
rz(1.0421853) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.790697) q[0];
sx q[0];
rz(-2.6804374) q[0];
sx q[0];
rz(2.8275369) q[0];
rz(0.028388609) q[1];
sx q[1];
rz(-0.24463618) q[1];
sx q[1];
rz(-1.8520741) q[1];
rz(-2.6311572) q[2];
sx q[2];
rz(-1.0344997) q[2];
sx q[2];
rz(-2.79881) q[2];
rz(2.3343347) q[3];
sx q[3];
rz(-1.1716598) q[3];
sx q[3];
rz(-1.4902478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
