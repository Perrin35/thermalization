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
rz(-1.550996) q[1];
sx q[1];
rz(-1.6928033) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4690404) q[0];
sx q[0];
rz(-1.7069478) q[0];
sx q[0];
rz(-0.37369136) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7291387) q[2];
sx q[2];
rz(-1.9082365) q[2];
sx q[2];
rz(-0.98357633) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6054666) q[1];
sx q[1];
rz(-1.8768633) q[1];
sx q[1];
rz(2.0293959) q[1];
x q[2];
rz(-1.3643614) q[3];
sx q[3];
rz(-0.41461333) q[3];
sx q[3];
rz(-1.9569091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68580565) q[2];
sx q[2];
rz(-1.1263584) q[2];
sx q[2];
rz(2.8796999) q[2];
rz(-0.68771466) q[3];
sx q[3];
rz(-2.3759638) q[3];
sx q[3];
rz(0.31622893) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19551109) q[0];
sx q[0];
rz(-1.8571778) q[0];
sx q[0];
rz(0.84445697) q[0];
rz(2.5020245) q[1];
sx q[1];
rz(-0.25194672) q[1];
sx q[1];
rz(-1.2676988) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8728392) q[0];
sx q[0];
rz(-1.6913101) q[0];
sx q[0];
rz(-2.2655649) q[0];
rz(-pi) q[1];
rz(-0.67633171) q[2];
sx q[2];
rz(-1.4478292) q[2];
sx q[2];
rz(2.7148397) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.87593791) q[1];
sx q[1];
rz(-0.75232435) q[1];
sx q[1];
rz(-3.1347514) q[1];
rz(2.5540939) q[3];
sx q[3];
rz(-2.2031257) q[3];
sx q[3];
rz(-1.2459696) q[3];
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
rz(-0.31288475) q[3];
sx q[3];
rz(-0.94066921) q[3];
sx q[3];
rz(2.5488034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6784356) q[0];
sx q[0];
rz(-2.1364697) q[0];
sx q[0];
rz(0.67331719) q[0];
rz(-2.8939269) q[1];
sx q[1];
rz(-1.0177871) q[1];
sx q[1];
rz(-2.8376875) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16557753) q[0];
sx q[0];
rz(-1.0212693) q[0];
sx q[0];
rz(-0.59161341) q[0];
rz(2.3788358) q[2];
sx q[2];
rz(-1.0029364) q[2];
sx q[2];
rz(3.0936846) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.60982882) q[1];
sx q[1];
rz(-1.7103238) q[1];
sx q[1];
rz(3.082485) q[1];
x q[2];
rz(2.9274793) q[3];
sx q[3];
rz(-1.9823834) q[3];
sx q[3];
rz(-2.0579091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5725382) q[2];
sx q[2];
rz(-1.7819541) q[2];
sx q[2];
rz(0.016810091) q[2];
rz(-2.861764) q[3];
sx q[3];
rz(-2.7357416) q[3];
sx q[3];
rz(-0.60389891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0458321) q[0];
sx q[0];
rz(-1.1149167) q[0];
sx q[0];
rz(-1.260585) q[0];
rz(0.49651217) q[1];
sx q[1];
rz(-1.9910201) q[1];
sx q[1];
rz(-1.6373434) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6034318) q[0];
sx q[0];
rz(-1.0350845) q[0];
sx q[0];
rz(-0.49749048) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7563633) q[2];
sx q[2];
rz(-2.5653815) q[2];
sx q[2];
rz(2.6389337) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1671686) q[1];
sx q[1];
rz(-1.4054055) q[1];
sx q[1];
rz(-1.0252293) q[1];
rz(2.639719) q[3];
sx q[3];
rz(-1.8326559) q[3];
sx q[3];
rz(3.1002432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.090195) q[2];
sx q[2];
rz(-0.53956705) q[2];
sx q[2];
rz(2.9810193) q[2];
rz(1.4422902) q[3];
sx q[3];
rz(-1.620404) q[3];
sx q[3];
rz(-3.1062533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2546688) q[0];
sx q[0];
rz(-1.2441607) q[0];
sx q[0];
rz(0.96920335) q[0];
rz(1.2508378) q[1];
sx q[1];
rz(-0.67986095) q[1];
sx q[1];
rz(-2.9108237) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7353775) q[0];
sx q[0];
rz(-3.1054524) q[0];
sx q[0];
rz(-1.4213495) q[0];
rz(-pi) q[1];
rz(-1.9491558) q[2];
sx q[2];
rz(-1.5659815) q[2];
sx q[2];
rz(0.7418599) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.240162) q[1];
sx q[1];
rz(-1.104875) q[1];
sx q[1];
rz(-0.45216143) q[1];
rz(-pi) q[2];
rz(-1.9649361) q[3];
sx q[3];
rz(-1.1232383) q[3];
sx q[3];
rz(2.7209865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1988924) q[2];
sx q[2];
rz(-2.5653699) q[2];
sx q[2];
rz(2.8885941) q[2];
rz(-1.3544072) q[3];
sx q[3];
rz(-1.1043045) q[3];
sx q[3];
rz(1.838292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.81786466) q[0];
sx q[0];
rz(-0.67960328) q[0];
sx q[0];
rz(-3.0815109) q[0];
rz(2.5443063) q[1];
sx q[1];
rz(-0.92514721) q[1];
sx q[1];
rz(-0.047860535) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56592076) q[0];
sx q[0];
rz(-1.4788155) q[0];
sx q[0];
rz(-3.1264114) q[0];
x q[1];
rz(-0.20345236) q[2];
sx q[2];
rz(-1.5849539) q[2];
sx q[2];
rz(-2.3490273) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7737433) q[1];
sx q[1];
rz(-1.7715104) q[1];
sx q[1];
rz(1.9108921) q[1];
x q[2];
rz(-1.490363) q[3];
sx q[3];
rz(-1.24986) q[3];
sx q[3];
rz(1.3453573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4045599) q[2];
sx q[2];
rz(-2.3768349) q[2];
sx q[2];
rz(2.9743189) q[2];
rz(0.078992756) q[3];
sx q[3];
rz(-1.9588214) q[3];
sx q[3];
rz(-2.3110068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9623742) q[0];
sx q[0];
rz(-3.0308864) q[0];
sx q[0];
rz(-3.0141444) q[0];
rz(0.92370644) q[1];
sx q[1];
rz(-2.5741003) q[1];
sx q[1];
rz(-2.4097402) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8475721) q[0];
sx q[0];
rz(-0.56041996) q[0];
sx q[0];
rz(-0.27604485) q[0];
rz(0.62035608) q[2];
sx q[2];
rz(-2.6292703) q[2];
sx q[2];
rz(-0.27119766) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3819363) q[1];
sx q[1];
rz(-2.3664673) q[1];
sx q[1];
rz(2.0630906) q[1];
x q[2];
rz(0.50639373) q[3];
sx q[3];
rz(-2.3558338) q[3];
sx q[3];
rz(-0.23791152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7544516) q[2];
sx q[2];
rz(-2.0595422) q[2];
sx q[2];
rz(-0.88073909) q[2];
rz(0.292101) q[3];
sx q[3];
rz(-2.5160774) q[3];
sx q[3];
rz(0.97091278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77325118) q[0];
sx q[0];
rz(-1.0057978) q[0];
sx q[0];
rz(2.3604426) q[0];
rz(1.4190405) q[1];
sx q[1];
rz(-0.96839372) q[1];
sx q[1];
rz(2.9636545) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3969361) q[0];
sx q[0];
rz(-1.2632003) q[0];
sx q[0];
rz(1.164695) q[0];
rz(-2.2841112) q[2];
sx q[2];
rz(-0.65543398) q[2];
sx q[2];
rz(1.0873742) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8101681) q[1];
sx q[1];
rz(-1.5123864) q[1];
sx q[1];
rz(-1.0845137) q[1];
rz(-1.7573886) q[3];
sx q[3];
rz(-1.3114138) q[3];
sx q[3];
rz(-2.4566539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6302744) q[2];
sx q[2];
rz(-1.4354458) q[2];
sx q[2];
rz(1.1159631) q[2];
rz(-0.80439862) q[3];
sx q[3];
rz(-0.65785995) q[3];
sx q[3];
rz(0.72949725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6804752) q[0];
sx q[0];
rz(-2.456433) q[0];
sx q[0];
rz(2.6035736) q[0];
rz(-2.8021725) q[1];
sx q[1];
rz(-1.598204) q[1];
sx q[1];
rz(3.0982049) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.246286) q[0];
sx q[0];
rz(-3.0247547) q[0];
sx q[0];
rz(2.1645592) q[0];
rz(1.0203982) q[2];
sx q[2];
rz(-1.4417037) q[2];
sx q[2];
rz(1.0132524) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0108569) q[1];
sx q[1];
rz(-0.82243516) q[1];
sx q[1];
rz(-1.6793516) q[1];
rz(-pi) q[2];
rz(-1.9566986) q[3];
sx q[3];
rz(-2.1348955) q[3];
sx q[3];
rz(-0.38340873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.033120774) q[2];
sx q[2];
rz(-2.7816009) q[2];
sx q[2];
rz(-1.6610422) q[2];
rz(0.2964274) q[3];
sx q[3];
rz(-2.0015643) q[3];
sx q[3];
rz(-0.75291434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.4794889) q[0];
sx q[0];
rz(-2.7487553) q[0];
sx q[0];
rz(1.9774849) q[0];
rz(2.8083943) q[1];
sx q[1];
rz(-1.4771799) q[1];
sx q[1];
rz(-2.3936757) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40953982) q[0];
sx q[0];
rz(-0.68258572) q[0];
sx q[0];
rz(-0.48614283) q[0];
rz(-pi) q[1];
rz(-2.203132) q[2];
sx q[2];
rz(-1.7424287) q[2];
sx q[2];
rz(0.29612637) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1031706) q[1];
sx q[1];
rz(-1.6254679) q[1];
sx q[1];
rz(1.5780539) q[1];
x q[2];
rz(0.87769784) q[3];
sx q[3];
rz(-1.6060281) q[3];
sx q[3];
rz(0.30634634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7738639) q[2];
sx q[2];
rz(-1.4089971) q[2];
sx q[2];
rz(1.9749953) q[2];
rz(1.599132) q[3];
sx q[3];
rz(-1.5365994) q[3];
sx q[3];
rz(-2.0994073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.790697) q[0];
sx q[0];
rz(-2.6804374) q[0];
sx q[0];
rz(2.8275369) q[0];
rz(3.113204) q[1];
sx q[1];
rz(-2.8969565) q[1];
sx q[1];
rz(1.2895186) q[1];
rz(0.88282925) q[2];
sx q[2];
rz(-0.7226823) q[2];
sx q[2];
rz(-0.48831627) q[2];
rz(-0.52846626) q[3];
sx q[3];
rz(-0.8800284) q[3];
sx q[3];
rz(2.8660091) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
