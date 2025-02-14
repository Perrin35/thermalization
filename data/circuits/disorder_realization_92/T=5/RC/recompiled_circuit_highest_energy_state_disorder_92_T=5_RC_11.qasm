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
rz(-0.31779274) q[0];
rz(-2.3843482) q[1];
sx q[1];
rz(4.7321893) q[1];
sx q[1];
rz(7.7319747) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4690404) q[0];
sx q[0];
rz(-1.4346449) q[0];
sx q[0];
rz(0.37369136) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4223675) q[2];
sx q[2];
rz(-2.7701391) q[2];
sx q[2];
rz(1.4329706) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6054666) q[1];
sx q[1];
rz(-1.2647293) q[1];
sx q[1];
rz(-1.1121967) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9775544) q[3];
sx q[3];
rz(-1.6534605) q[3];
sx q[3];
rz(0.57549046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.68580565) q[2];
sx q[2];
rz(-2.0152342) q[2];
sx q[2];
rz(-2.8796999) q[2];
rz(0.68771466) q[3];
sx q[3];
rz(-0.7656289) q[3];
sx q[3];
rz(0.31622893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9460816) q[0];
sx q[0];
rz(-1.2844149) q[0];
sx q[0];
rz(2.2971357) q[0];
rz(2.5020245) q[1];
sx q[1];
rz(-2.8896459) q[1];
sx q[1];
rz(1.2676988) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2687534) q[0];
sx q[0];
rz(-1.6913101) q[0];
sx q[0];
rz(0.87602776) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67633171) q[2];
sx q[2];
rz(-1.4478292) q[2];
sx q[2];
rz(2.7148397) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2750249) q[1];
sx q[1];
rz(-2.3230988) q[1];
sx q[1];
rz(-1.5643934) q[1];
rz(2.5540939) q[3];
sx q[3];
rz(-2.2031257) q[3];
sx q[3];
rz(1.895623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1203479) q[2];
sx q[2];
rz(-2.6159888) q[2];
sx q[2];
rz(2.8525103) q[2];
rz(-0.31288475) q[3];
sx q[3];
rz(-0.94066921) q[3];
sx q[3];
rz(2.5488034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6784356) q[0];
sx q[0];
rz(-1.005123) q[0];
sx q[0];
rz(2.4682755) q[0];
rz(-0.24766573) q[1];
sx q[1];
rz(-1.0177871) q[1];
sx q[1];
rz(2.8376875) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0659828) q[0];
sx q[0];
rz(-2.3572266) q[0];
sx q[0];
rz(0.83215587) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2938305) q[2];
sx q[2];
rz(-0.94901949) q[2];
sx q[2];
rz(-1.997681) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.188855) q[1];
sx q[1];
rz(-1.6293289) q[1];
sx q[1];
rz(-1.7105647) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.21411333) q[3];
sx q[3];
rz(-1.1592093) q[3];
sx q[3];
rz(2.0579091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5725382) q[2];
sx q[2];
rz(-1.7819541) q[2];
sx q[2];
rz(-3.1247826) q[2];
rz(2.861764) q[3];
sx q[3];
rz(-0.40585104) q[3];
sx q[3];
rz(-0.60389891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0957606) q[0];
sx q[0];
rz(-1.1149167) q[0];
sx q[0];
rz(1.260585) q[0];
rz(-0.49651217) q[1];
sx q[1];
rz(-1.1505726) q[1];
sx q[1];
rz(1.5042492) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8634254) q[0];
sx q[0];
rz(-2.4274732) q[0];
sx q[0];
rz(-0.89366736) q[0];
x q[1];
rz(1.331318) q[2];
sx q[2];
rz(-1.0415072) q[2];
sx q[2];
rz(-3.0893507) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.28037) q[1];
sx q[1];
rz(-2.5739453) q[1];
sx q[1];
rz(-1.2595792) q[1];
rz(-1.8674866) q[3];
sx q[3];
rz(-1.0875351) q[3];
sx q[3];
rz(-1.7532574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.090195) q[2];
sx q[2];
rz(-0.53956705) q[2];
sx q[2];
rz(0.16057333) q[2];
rz(1.4422902) q[3];
sx q[3];
rz(-1.5211886) q[3];
sx q[3];
rz(-0.035339385) q[3];
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
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8869239) q[0];
sx q[0];
rz(-1.897432) q[0];
sx q[0];
rz(0.96920335) q[0];
rz(-1.8907549) q[1];
sx q[1];
rz(-0.67986095) q[1];
sx q[1];
rz(0.23076898) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8849205) q[0];
sx q[0];
rz(-1.6065336) q[0];
sx q[0];
rz(3.1362094) q[0];
rz(1.5838301) q[2];
sx q[2];
rz(-0.37838867) q[2];
sx q[2];
rz(-2.3005444) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0645551) q[1];
sx q[1];
rz(-2.5041576) q[1];
sx q[1];
rz(2.2861479) q[1];
x q[2];
rz(1.1766566) q[3];
sx q[3];
rz(-1.1232383) q[3];
sx q[3];
rz(2.7209865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.9427003) q[2];
sx q[2];
rz(-0.57622272) q[2];
sx q[2];
rz(-2.8885941) q[2];
rz(1.7871855) q[3];
sx q[3];
rz(-2.0372882) q[3];
sx q[3];
rz(-1.838292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.323728) q[0];
sx q[0];
rz(-2.4619894) q[0];
sx q[0];
rz(-0.060081765) q[0];
rz(-0.59728638) q[1];
sx q[1];
rz(-2.2164454) q[1];
sx q[1];
rz(0.047860535) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0062701) q[0];
sx q[0];
rz(-1.5859134) q[0];
sx q[0];
rz(-1.478805) q[0];
rz(-pi) q[1];
x q[1];
rz(1.585252) q[2];
sx q[2];
rz(-1.3673646) q[2];
sx q[2];
rz(-2.360441) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4255449) q[1];
sx q[1];
rz(-2.7486781) q[1];
sx q[1];
rz(-1.0231189) q[1];
rz(-pi) q[2];
rz(1.490363) q[3];
sx q[3];
rz(-1.24986) q[3];
sx q[3];
rz(-1.3453573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.73703274) q[2];
sx q[2];
rz(-0.76475778) q[2];
sx q[2];
rz(-0.16727373) q[2];
rz(-0.078992756) q[3];
sx q[3];
rz(-1.1827712) q[3];
sx q[3];
rz(-2.3110068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9623742) q[0];
sx q[0];
rz(-0.11070624) q[0];
sx q[0];
rz(0.12744823) q[0];
rz(-0.92370644) q[1];
sx q[1];
rz(-2.5741003) q[1];
sx q[1];
rz(-0.73185241) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028721044) q[0];
sx q[0];
rz(-2.1076308) q[0];
sx q[0];
rz(-1.4013994) q[0];
rz(-1.8867854) q[2];
sx q[2];
rz(-1.1605186) q[2];
sx q[2];
rz(-0.41543322) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7596563) q[1];
sx q[1];
rz(-0.77512533) q[1];
sx q[1];
rz(1.0785021) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1189156) q[3];
sx q[3];
rz(-2.2377399) q[3];
sx q[3];
rz(-0.42740145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.38714108) q[2];
sx q[2];
rz(-1.0820505) q[2];
sx q[2];
rz(0.88073909) q[2];
rz(2.8494917) q[3];
sx q[3];
rz(-2.5160774) q[3];
sx q[3];
rz(2.1706799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77325118) q[0];
sx q[0];
rz(-1.0057978) q[0];
sx q[0];
rz(-0.7811501) q[0];
rz(-1.4190405) q[1];
sx q[1];
rz(-0.96839372) q[1];
sx q[1];
rz(-2.9636545) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74465655) q[0];
sx q[0];
rz(-1.8783924) q[0];
sx q[0];
rz(-1.9768977) q[0];
rz(-pi) q[1];
rz(2.2841112) q[2];
sx q[2];
rz(-0.65543398) q[2];
sx q[2];
rz(-1.0873742) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3314245) q[1];
sx q[1];
rz(-1.5123864) q[1];
sx q[1];
rz(2.057079) q[1];
rz(0.61011647) q[3];
sx q[3];
rz(-0.31829208) q[3];
sx q[3];
rz(1.3194609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6302744) q[2];
sx q[2];
rz(-1.7061468) q[2];
sx q[2];
rz(-2.0256296) q[2];
rz(2.337194) q[3];
sx q[3];
rz(-2.4837327) q[3];
sx q[3];
rz(-0.72949725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6804752) q[0];
sx q[0];
rz(-0.68515968) q[0];
sx q[0];
rz(-0.53801909) q[0];
rz(-2.8021725) q[1];
sx q[1];
rz(-1.5433886) q[1];
sx q[1];
rz(-3.0982049) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8432253) q[0];
sx q[0];
rz(-1.6675673) q[0];
sx q[0];
rz(-3.0760187) q[0];
rz(2.9904463) q[2];
sx q[2];
rz(-1.0254964) q[2];
sx q[2];
rz(0.4787094) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.85203899) q[1];
sx q[1];
rz(-2.3869136) q[1];
sx q[1];
rz(-0.11615495) q[1];
rz(2.6049117) q[3];
sx q[3];
rz(-2.4701475) q[3];
sx q[3];
rz(1.0332024) q[3];
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
rz(-1.4805504) q[2];
rz(-0.2964274) q[3];
sx q[3];
rz(-2.0015643) q[3];
sx q[3];
rz(0.75291434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66210371) q[0];
sx q[0];
rz(-2.7487553) q[0];
sx q[0];
rz(-1.9774849) q[0];
rz(0.33319831) q[1];
sx q[1];
rz(-1.6644128) q[1];
sx q[1];
rz(-2.3936757) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5503905) q[0];
sx q[0];
rz(-1.2716312) q[0];
sx q[0];
rz(-2.5183866) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9299253) q[2];
sx q[2];
rz(-0.94918409) q[2];
sx q[2];
rz(1.399161) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.038422) q[1];
sx q[1];
rz(-1.5161247) q[1];
sx q[1];
rz(-1.5780539) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2638948) q[3];
sx q[3];
rz(-1.5355645) q[3];
sx q[3];
rz(-0.30634634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7738639) q[2];
sx q[2];
rz(-1.7325956) q[2];
sx q[2];
rz(-1.9749953) q[2];
rz(-1.599132) q[3];
sx q[3];
rz(-1.5365994) q[3];
sx q[3];
rz(-1.0421853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3508956) q[0];
sx q[0];
rz(-2.6804374) q[0];
sx q[0];
rz(2.8275369) q[0];
rz(-0.028388609) q[1];
sx q[1];
rz(-2.8969565) q[1];
sx q[1];
rz(1.2895186) q[1];
rz(-2.1688228) q[2];
sx q[2];
rz(-1.1373873) q[2];
sx q[2];
rz(1.6349229) q[2];
rz(-2.1185067) q[3];
sx q[3];
rz(-2.2991091) q[3];
sx q[3];
rz(-2.6753796) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
