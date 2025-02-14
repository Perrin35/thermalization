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
rz(-2.3075624) q[0];
sx q[0];
rz(-0.40606719) q[0];
sx q[0];
rz(1.1817482) q[0];
rz(1.740068) q[1];
sx q[1];
rz(5.8086173) q[1];
sx q[1];
rz(9.557815) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1057464) q[0];
sx q[0];
rz(-1.3091358) q[0];
sx q[0];
rz(-0.91069551) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70212097) q[2];
sx q[2];
rz(-0.63533855) q[2];
sx q[2];
rz(1.8566657) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0892355) q[1];
sx q[1];
rz(-1.3468139) q[1];
sx q[1];
rz(1.2473945) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0361133) q[3];
sx q[3];
rz(-0.64247349) q[3];
sx q[3];
rz(1.5657263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0066321) q[2];
sx q[2];
rz(-2.1217608) q[2];
sx q[2];
rz(0.68438619) q[2];
rz(2.0002174) q[3];
sx q[3];
rz(-0.26641521) q[3];
sx q[3];
rz(-0.086708955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.3619277) q[0];
sx q[0];
rz(-0.68988887) q[0];
sx q[0];
rz(-0.46098125) q[0];
rz(1.1191248) q[1];
sx q[1];
rz(-2.7179167) q[1];
sx q[1];
rz(-0.31203312) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082816007) q[0];
sx q[0];
rz(-0.030555336) q[0];
sx q[0];
rz(-2.729134) q[0];
rz(0.076032555) q[2];
sx q[2];
rz(-0.8113779) q[2];
sx q[2];
rz(-1.8234314) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4496135) q[1];
sx q[1];
rz(-0.17573729) q[1];
sx q[1];
rz(-1.537561) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8617763) q[3];
sx q[3];
rz(-0.93799611) q[3];
sx q[3];
rz(2.1386843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.95388874) q[2];
sx q[2];
rz(-1.1053332) q[2];
sx q[2];
rz(0.40404955) q[2];
rz(-0.36014253) q[3];
sx q[3];
rz(-1.4934243) q[3];
sx q[3];
rz(-0.68296105) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78524154) q[0];
sx q[0];
rz(-1.0026362) q[0];
sx q[0];
rz(0.32401618) q[0];
rz(2.0815381) q[1];
sx q[1];
rz(-2.8476604) q[1];
sx q[1];
rz(-0.70045984) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9793036) q[0];
sx q[0];
rz(-1.7320181) q[0];
sx q[0];
rz(0.17031807) q[0];
rz(-pi) q[1];
rz(1.9159928) q[2];
sx q[2];
rz(-1.7734062) q[2];
sx q[2];
rz(2.4576296) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6571275) q[1];
sx q[1];
rz(-1.8052237) q[1];
sx q[1];
rz(2.3695494) q[1];
x q[2];
rz(1.1491997) q[3];
sx q[3];
rz(-1.7497049) q[3];
sx q[3];
rz(-0.94551555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5141653) q[2];
sx q[2];
rz(-2.8943987) q[2];
sx q[2];
rz(-2.3449281) q[2];
rz(1.1967777) q[3];
sx q[3];
rz(-1.6007042) q[3];
sx q[3];
rz(0.59320199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6881123) q[0];
sx q[0];
rz(-1.0013094) q[0];
sx q[0];
rz(1.8185115) q[0];
rz(0.46609136) q[1];
sx q[1];
rz(-0.90704647) q[1];
sx q[1];
rz(-1.1493433) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2043484) q[0];
sx q[0];
rz(-0.68182985) q[0];
sx q[0];
rz(2.3005062) q[0];
rz(-0.39300275) q[2];
sx q[2];
rz(-2.1519024) q[2];
sx q[2];
rz(-0.71497922) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6862029) q[1];
sx q[1];
rz(-1.7727994) q[1];
sx q[1];
rz(-1.6167902) q[1];
x q[2];
rz(-2.9907865) q[3];
sx q[3];
rz(-1.9016087) q[3];
sx q[3];
rz(2.2589333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.65679067) q[2];
sx q[2];
rz(-0.66978407) q[2];
sx q[2];
rz(2.316324) q[2];
rz(-1.2800062) q[3];
sx q[3];
rz(-1.5214336) q[3];
sx q[3];
rz(-2.4204204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.698302) q[0];
sx q[0];
rz(-1.8809141) q[0];
sx q[0];
rz(1.3662421) q[0];
rz(-1.0900991) q[1];
sx q[1];
rz(-1.6366199) q[1];
sx q[1];
rz(-1.3390138) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1813359) q[0];
sx q[0];
rz(-2.4614545) q[0];
sx q[0];
rz(-2.8614317) q[0];
rz(-0.86651643) q[2];
sx q[2];
rz(-2.0712059) q[2];
sx q[2];
rz(1.6472218) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.80833188) q[1];
sx q[1];
rz(-1.3587316) q[1];
sx q[1];
rz(1.6912837) q[1];
rz(0.7991661) q[3];
sx q[3];
rz(-0.46119565) q[3];
sx q[3];
rz(0.14652625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.94685405) q[2];
sx q[2];
rz(-1.4118492) q[2];
sx q[2];
rz(-0.33352387) q[2];
rz(0.58878318) q[3];
sx q[3];
rz(-2.6684561) q[3];
sx q[3];
rz(0.80140448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81353417) q[0];
sx q[0];
rz(-1.3575587) q[0];
sx q[0];
rz(1.300977) q[0];
rz(-2.1133555) q[1];
sx q[1];
rz(-2.6515549) q[1];
sx q[1];
rz(-2.6992544) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16602941) q[0];
sx q[0];
rz(-1.1021758) q[0];
sx q[0];
rz(1.2965078) q[0];
rz(2.2919334) q[2];
sx q[2];
rz(-2.4772212) q[2];
sx q[2];
rz(-2.0400782) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.42980151) q[1];
sx q[1];
rz(-0.45908976) q[1];
sx q[1];
rz(2.7613822) q[1];
rz(-0.30486561) q[3];
sx q[3];
rz(-2.8072522) q[3];
sx q[3];
rz(-1.6337194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3698795) q[2];
sx q[2];
rz(-1.5481202) q[2];
sx q[2];
rz(-2.6698574) q[2];
rz(-1.6996023) q[3];
sx q[3];
rz(-2.580018) q[3];
sx q[3];
rz(2.877318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0278552) q[0];
sx q[0];
rz(-1.6738482) q[0];
sx q[0];
rz(-0.16673985) q[0];
rz(-0.79404229) q[1];
sx q[1];
rz(-1.200095) q[1];
sx q[1];
rz(-3.0418495) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82631153) q[0];
sx q[0];
rz(-2.4874685) q[0];
sx q[0];
rz(-2.099206) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0641618) q[2];
sx q[2];
rz(-0.68175137) q[2];
sx q[2];
rz(-0.1296986) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.9445716) q[1];
sx q[1];
rz(-0.58496956) q[1];
sx q[1];
rz(2.9119063) q[1];
rz(-pi) q[2];
rz(2.0896627) q[3];
sx q[3];
rz(-0.73118932) q[3];
sx q[3];
rz(2.2377781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2793067) q[2];
sx q[2];
rz(-2.0760832) q[2];
sx q[2];
rz(3.0118937) q[2];
rz(1.2636403) q[3];
sx q[3];
rz(-1.940515) q[3];
sx q[3];
rz(-2.9960846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97067863) q[0];
sx q[0];
rz(-2.2190974) q[0];
sx q[0];
rz(-2.7150735) q[0];
rz(1.0921987) q[1];
sx q[1];
rz(-2.009232) q[1];
sx q[1];
rz(-2.138864) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20518453) q[0];
sx q[0];
rz(-0.6615839) q[0];
sx q[0];
rz(-1.9993881) q[0];
x q[1];
rz(-0.10522225) q[2];
sx q[2];
rz(-2.4159523) q[2];
sx q[2];
rz(1.2187472) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7510609) q[1];
sx q[1];
rz(-0.56274429) q[1];
sx q[1];
rz(-2.2319002) q[1];
rz(-pi) q[2];
rz(1.1158607) q[3];
sx q[3];
rz(-1.7835894) q[3];
sx q[3];
rz(1.7011332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2670474) q[2];
sx q[2];
rz(-2.5858877) q[2];
sx q[2];
rz(0.12346539) q[2];
rz(-0.63001436) q[3];
sx q[3];
rz(-1.29653) q[3];
sx q[3];
rz(0.71648487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89710871) q[0];
sx q[0];
rz(-0.69863313) q[0];
sx q[0];
rz(0.83936349) q[0];
rz(-2.9821303) q[1];
sx q[1];
rz(-2.1493252) q[1];
sx q[1];
rz(-2.2119567) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2004505) q[0];
sx q[0];
rz(-1.3172035) q[0];
sx q[0];
rz(-2.3932049) q[0];
rz(2.9603297) q[2];
sx q[2];
rz(-0.59872228) q[2];
sx q[2];
rz(1.8229654) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8690455) q[1];
sx q[1];
rz(-2.2490977) q[1];
sx q[1];
rz(2.9806304) q[1];
rz(0.93299039) q[3];
sx q[3];
rz(-1.3853042) q[3];
sx q[3];
rz(0.18292566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.9440426) q[2];
sx q[2];
rz(-1.1669179) q[2];
sx q[2];
rz(2.3738532) q[2];
rz(-2.7770216) q[3];
sx q[3];
rz(-1.3837827) q[3];
sx q[3];
rz(0.069469623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1345054) q[0];
sx q[0];
rz(-2.5275079) q[0];
sx q[0];
rz(-2.2579204) q[0];
rz(2.9332352) q[1];
sx q[1];
rz(-0.50269428) q[1];
sx q[1];
rz(1.8066033) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3985838) q[0];
sx q[0];
rz(-1.5129876) q[0];
sx q[0];
rz(3.0763118) q[0];
rz(1.7572549) q[2];
sx q[2];
rz(-1.0508732) q[2];
sx q[2];
rz(-2.7911719) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61858863) q[1];
sx q[1];
rz(-0.72806137) q[1];
sx q[1];
rz(0.87597202) q[1];
rz(-pi) q[2];
rz(2.317071) q[3];
sx q[3];
rz(-2.093165) q[3];
sx q[3];
rz(2.8461412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9779382) q[2];
sx q[2];
rz(-1.5074707) q[2];
sx q[2];
rz(3.0955691) q[2];
rz(2.9764791) q[3];
sx q[3];
rz(-2.8486227) q[3];
sx q[3];
rz(-1.2142115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5788427) q[0];
sx q[0];
rz(-3.0414707) q[0];
sx q[0];
rz(1.1895251) q[0];
rz(2.9764755) q[1];
sx q[1];
rz(-1.4585635) q[1];
sx q[1];
rz(-0.30452902) q[1];
rz(2.5223059) q[2];
sx q[2];
rz(-2.613667) q[2];
sx q[2];
rz(2.1783301) q[2];
rz(-0.48043171) q[3];
sx q[3];
rz(-1.1050925) q[3];
sx q[3];
rz(2.1375124) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
