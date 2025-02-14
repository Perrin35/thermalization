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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6725522) q[0];
sx q[0];
rz(-1.7069478) q[0];
sx q[0];
rz(-0.37369136) q[0];
x q[1];
rz(0.4223675) q[2];
sx q[2];
rz(-2.7701391) q[2];
sx q[2];
rz(1.708622) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9592465) q[1];
sx q[1];
rz(-1.1350147) q[1];
sx q[1];
rz(-2.8027727) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1640383) q[3];
sx q[3];
rz(-1.4881322) q[3];
sx q[3];
rz(2.5661022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.68580565) q[2];
sx q[2];
rz(-1.1263584) q[2];
sx q[2];
rz(-0.2618928) q[2];
rz(-0.68771466) q[3];
sx q[3];
rz(-2.3759638) q[3];
sx q[3];
rz(0.31622893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9460816) q[0];
sx q[0];
rz(-1.8571778) q[0];
sx q[0];
rz(-0.84445697) q[0];
rz(0.63956815) q[1];
sx q[1];
rz(-0.25194672) q[1];
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
rz(0.19494178) q[2];
sx q[2];
rz(-2.4559074) q[2];
sx q[2];
rz(0.99239381) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.86656777) q[1];
sx q[1];
rz(-0.81849388) q[1];
sx q[1];
rz(1.5771992) q[1];
x q[2];
rz(0.92314641) q[3];
sx q[3];
rz(-0.83469838) q[3];
sx q[3];
rz(-2.0908963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1203479) q[2];
sx q[2];
rz(-0.52560386) q[2];
sx q[2];
rz(-0.28908238) q[2];
rz(2.8287079) q[3];
sx q[3];
rz(-0.94066921) q[3];
sx q[3];
rz(2.5488034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4631571) q[0];
sx q[0];
rz(-1.005123) q[0];
sx q[0];
rz(2.4682755) q[0];
rz(-2.8939269) q[1];
sx q[1];
rz(-1.0177871) q[1];
sx q[1];
rz(0.30390513) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.073868) q[0];
sx q[0];
rz(-1.0751599) q[0];
sx q[0];
rz(0.93511029) q[0];
x q[1];
rz(-2.2938305) q[2];
sx q[2];
rz(-0.94901949) q[2];
sx q[2];
rz(-1.997681) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60982882) q[1];
sx q[1];
rz(-1.7103238) q[1];
sx q[1];
rz(3.082485) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1506733) q[3];
sx q[3];
rz(-1.3748079) q[3];
sx q[3];
rz(0.57388692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5725382) q[2];
sx q[2];
rz(-1.3596386) q[2];
sx q[2];
rz(-0.016810091) q[2];
rz(-0.2798287) q[3];
sx q[3];
rz(-2.7357416) q[3];
sx q[3];
rz(0.60389891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.0458321) q[0];
sx q[0];
rz(-2.026676) q[0];
sx q[0];
rz(1.260585) q[0];
rz(-2.6450805) q[1];
sx q[1];
rz(-1.1505726) q[1];
sx q[1];
rz(-1.5042492) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53816089) q[0];
sx q[0];
rz(-2.1065082) q[0];
sx q[0];
rz(-2.6441022) q[0];
x q[1];
rz(-2.7563633) q[2];
sx q[2];
rz(-2.5653815) q[2];
sx q[2];
rz(2.6389337) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.97442409) q[1];
sx q[1];
rz(-1.4054055) q[1];
sx q[1];
rz(-2.1163634) q[1];
x q[2];
rz(0.50187364) q[3];
sx q[3];
rz(-1.3089367) q[3];
sx q[3];
rz(3.1002432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.090195) q[2];
sx q[2];
rz(-2.6020256) q[2];
sx q[2];
rz(2.9810193) q[2];
rz(-1.6993025) q[3];
sx q[3];
rz(-1.620404) q[3];
sx q[3];
rz(-3.1062533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2546688) q[0];
sx q[0];
rz(-1.2441607) q[0];
sx q[0];
rz(-2.1723893) q[0];
rz(-1.2508378) q[1];
sx q[1];
rz(-0.67986095) q[1];
sx q[1];
rz(2.9108237) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7353775) q[0];
sx q[0];
rz(-3.1054524) q[0];
sx q[0];
rz(1.4213495) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0051813263) q[2];
sx q[2];
rz(-1.1924414) q[2];
sx q[2];
rz(0.82702247) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2573852) q[1];
sx q[1];
rz(-1.9717968) q[1];
sx q[1];
rz(2.0805418) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9649361) q[3];
sx q[3];
rz(-2.0183544) q[3];
sx q[3];
rz(0.42060619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1988924) q[2];
sx q[2];
rz(-0.57622272) q[2];
sx q[2];
rz(-0.25299859) q[2];
rz(-1.7871855) q[3];
sx q[3];
rz(-2.0372882) q[3];
sx q[3];
rz(1.838292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.323728) q[0];
sx q[0];
rz(-2.4619894) q[0];
sx q[0];
rz(-3.0815109) q[0];
rz(0.59728638) q[1];
sx q[1];
rz(-2.2164454) q[1];
sx q[1];
rz(-0.047860535) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7394848) q[0];
sx q[0];
rz(-3.0483709) q[0];
sx q[0];
rz(-1.4076821) q[0];
rz(-pi) q[1];
rz(1.5563407) q[2];
sx q[2];
rz(-1.774228) q[2];
sx q[2];
rz(0.78115168) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.3678494) q[1];
sx q[1];
rz(-1.7715104) q[1];
sx q[1];
rz(-1.2307005) q[1];
rz(0.32190692) q[3];
sx q[3];
rz(-1.4944781) q[3];
sx q[3];
rz(-2.8907311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73703274) q[2];
sx q[2];
rz(-0.76475778) q[2];
sx q[2];
rz(2.9743189) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(1.5122685) q[0];
sx q[0];
rz(-1.716181) q[0];
sx q[0];
rz(0.54319197) q[0];
rz(-0.62035608) q[2];
sx q[2];
rz(-2.6292703) q[2];
sx q[2];
rz(0.27119766) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5868581) q[1];
sx q[1];
rz(-1.233685) q[1];
sx q[1];
rz(2.2829775) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1189156) q[3];
sx q[3];
rz(-2.2377399) q[3];
sx q[3];
rz(-2.7141912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7544516) q[2];
sx q[2];
rz(-2.0595422) q[2];
sx q[2];
rz(0.88073909) q[2];
rz(0.292101) q[3];
sx q[3];
rz(-0.62551522) q[3];
sx q[3];
rz(2.1706799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77325118) q[0];
sx q[0];
rz(-1.0057978) q[0];
sx q[0];
rz(0.7811501) q[0];
rz(1.7225522) q[1];
sx q[1];
rz(-0.96839372) q[1];
sx q[1];
rz(-2.9636545) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3969361) q[0];
sx q[0];
rz(-1.2632003) q[0];
sx q[0];
rz(-1.9768977) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0441842) q[2];
sx q[2];
rz(-1.9810314) q[2];
sx q[2];
rz(-3.0237449) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8101681) q[1];
sx q[1];
rz(-1.5123864) q[1];
sx q[1];
rz(2.057079) q[1];
rz(-pi) q[2];
x q[2];
rz(1.384204) q[3];
sx q[3];
rz(-1.3114138) q[3];
sx q[3];
rz(0.6849388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.51131821) q[2];
sx q[2];
rz(-1.4354458) q[2];
sx q[2];
rz(2.0256296) q[2];
rz(2.337194) q[3];
sx q[3];
rz(-0.65785995) q[3];
sx q[3];
rz(0.72949725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46111742) q[0];
sx q[0];
rz(-0.68515968) q[0];
sx q[0];
rz(-0.53801909) q[0];
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
rz(0.26608426) q[0];
sx q[0];
rz(-1.5055297) q[0];
sx q[0];
rz(-1.6677744) q[0];
rz(-0.15114637) q[2];
sx q[2];
rz(-1.0254964) q[2];
sx q[2];
rz(0.4787094) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0108569) q[1];
sx q[1];
rz(-0.82243516) q[1];
sx q[1];
rz(1.6793516) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6049117) q[3];
sx q[3];
rz(-2.4701475) q[3];
sx q[3];
rz(1.0332024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.033120774) q[2];
sx q[2];
rz(-2.7816009) q[2];
sx q[2];
rz(1.6610422) q[2];
rz(-0.2964274) q[3];
sx q[3];
rz(-2.0015643) q[3];
sx q[3];
rz(-2.3886783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66210371) q[0];
sx q[0];
rz(-2.7487553) q[0];
sx q[0];
rz(-1.1641077) q[0];
rz(-0.33319831) q[1];
sx q[1];
rz(-1.4771799) q[1];
sx q[1];
rz(-2.3936757) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7320528) q[0];
sx q[0];
rz(-0.68258572) q[0];
sx q[0];
rz(2.6554498) q[0];
rz(-pi) q[1];
rz(1.8560772) q[2];
sx q[2];
rz(-0.65212265) q[2];
sx q[2];
rz(-1.0456606) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.235215) q[1];
sx q[1];
rz(-3.0864419) q[1];
sx q[1];
rz(-0.13184594) q[1];
rz(-0.045785964) q[3];
sx q[3];
rz(-0.87821315) q[3];
sx q[3];
rz(1.847895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.36772874) q[2];
sx q[2];
rz(-1.7325956) q[2];
sx q[2];
rz(1.9749953) q[2];
rz(-1.5424607) q[3];
sx q[3];
rz(-1.5365994) q[3];
sx q[3];
rz(1.0421853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3508956) q[0];
sx q[0];
rz(-0.46115524) q[0];
sx q[0];
rz(-0.3140558) q[0];
rz(3.113204) q[1];
sx q[1];
rz(-2.8969565) q[1];
sx q[1];
rz(1.2895186) q[1];
rz(0.97276982) q[2];
sx q[2];
rz(-1.1373873) q[2];
sx q[2];
rz(1.6349229) q[2];
rz(-2.3343347) q[3];
sx q[3];
rz(-1.9699329) q[3];
sx q[3];
rz(1.6513449) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
