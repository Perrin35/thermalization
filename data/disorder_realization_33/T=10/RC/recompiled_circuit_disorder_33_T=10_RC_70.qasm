OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5053951) q[0];
sx q[0];
rz(-2.8656821) q[0];
sx q[0];
rz(-1.3077868) q[0];
rz(1.1360599) q[1];
sx q[1];
rz(-0.93568957) q[1];
sx q[1];
rz(-1.5712665) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64496541) q[0];
sx q[0];
rz(-2.6500406) q[0];
sx q[0];
rz(0.77792032) q[0];
rz(-pi) q[1];
rz(0.88408006) q[2];
sx q[2];
rz(-0.93187098) q[2];
sx q[2];
rz(-2.0761348) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.624144) q[1];
sx q[1];
rz(-1.9074719) q[1];
sx q[1];
rz(1.3144073) q[1];
rz(-pi) q[2];
rz(-1.2075495) q[3];
sx q[3];
rz(-1.3875291) q[3];
sx q[3];
rz(-2.7024384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87542614) q[2];
sx q[2];
rz(-0.29310075) q[2];
sx q[2];
rz(2.0092633) q[2];
rz(1.4663565) q[3];
sx q[3];
rz(-1.8050067) q[3];
sx q[3];
rz(-2.1291389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9448626) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(-0.18584132) q[0];
rz(-0.56022412) q[1];
sx q[1];
rz(-1.2954243) q[1];
sx q[1];
rz(-2.9247608) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8502055) q[0];
sx q[0];
rz(-2.4225525) q[0];
sx q[0];
rz(2.015381) q[0];
rz(-pi) q[1];
rz(-2.6685643) q[2];
sx q[2];
rz(-1.066726) q[2];
sx q[2];
rz(2.8216528) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20570457) q[1];
sx q[1];
rz(-1.7824031) q[1];
sx q[1];
rz(0.88797027) q[1];
rz(-pi) q[2];
rz(0.84083765) q[3];
sx q[3];
rz(-2.3008122) q[3];
sx q[3];
rz(-0.18883146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8314787) q[2];
sx q[2];
rz(-2.3159413) q[2];
sx q[2];
rz(-1.8537834) q[2];
rz(-0.76256049) q[3];
sx q[3];
rz(-1.972714) q[3];
sx q[3];
rz(-2.8365703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6771616) q[0];
sx q[0];
rz(-2.7966249) q[0];
sx q[0];
rz(2.537354) q[0];
rz(-1.3263946) q[1];
sx q[1];
rz(-1.3605958) q[1];
sx q[1];
rz(2.2089829) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37768294) q[0];
sx q[0];
rz(-1.8881067) q[0];
sx q[0];
rz(-3.0850287) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2432125) q[2];
sx q[2];
rz(-0.2873688) q[2];
sx q[2];
rz(0.76295602) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.052913594) q[1];
sx q[1];
rz(-1.0148078) q[1];
sx q[1];
rz(1.4312137) q[1];
x q[2];
rz(-1.014939) q[3];
sx q[3];
rz(-1.1851289) q[3];
sx q[3];
rz(-0.25394299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9937667) q[2];
sx q[2];
rz(-2.0596762) q[2];
sx q[2];
rz(-1.0926584) q[2];
rz(-0.5422194) q[3];
sx q[3];
rz(-2.0565624) q[3];
sx q[3];
rz(-2.1742163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3595235) q[0];
sx q[0];
rz(-0.096465915) q[0];
sx q[0];
rz(0.50022593) q[0];
rz(-0.80530986) q[1];
sx q[1];
rz(-1.1601245) q[1];
sx q[1];
rz(-1.6436228) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1456137) q[0];
sx q[0];
rz(-2.0659475) q[0];
sx q[0];
rz(0.33546319) q[0];
rz(-pi) q[1];
rz(-0.056604071) q[2];
sx q[2];
rz(-1.7610234) q[2];
sx q[2];
rz(2.9375926) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8151617) q[1];
sx q[1];
rz(-1.8306499) q[1];
sx q[1];
rz(1.3304779) q[1];
rz(2.3718932) q[3];
sx q[3];
rz(-0.61004988) q[3];
sx q[3];
rz(-1.1806012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.74636373) q[2];
sx q[2];
rz(-0.56240288) q[2];
sx q[2];
rz(2.4397819) q[2];
rz(0.83135215) q[3];
sx q[3];
rz(-0.96389198) q[3];
sx q[3];
rz(-2.519616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2410626) q[0];
sx q[0];
rz(-2.5456972) q[0];
sx q[0];
rz(2.3262614) q[0];
rz(-1.6197846) q[1];
sx q[1];
rz(-0.83414572) q[1];
sx q[1];
rz(-2.0933847) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76083714) q[0];
sx q[0];
rz(-1.8341944) q[0];
sx q[0];
rz(0.25816985) q[0];
x q[1];
rz(-1.5830718) q[2];
sx q[2];
rz(-0.91931146) q[2];
sx q[2];
rz(-2.2001681) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.293562) q[1];
sx q[1];
rz(-1.8082431) q[1];
sx q[1];
rz(0.32229396) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1225554) q[3];
sx q[3];
rz(-0.69283797) q[3];
sx q[3];
rz(0.90941959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.52577019) q[2];
sx q[2];
rz(-2.5746391) q[2];
sx q[2];
rz(1.099951) q[2];
rz(-2.3163017) q[3];
sx q[3];
rz(-2.0992978) q[3];
sx q[3];
rz(2.2560789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0734171) q[0];
sx q[0];
rz(-0.59403479) q[0];
sx q[0];
rz(2.2391879) q[0];
rz(-1.0166608) q[1];
sx q[1];
rz(-1.0598176) q[1];
sx q[1];
rz(0.12983233) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8311365) q[0];
sx q[0];
rz(-1.2015011) q[0];
sx q[0];
rz(-2.5184758) q[0];
rz(-pi) q[1];
rz(2.1075222) q[2];
sx q[2];
rz(-2.495129) q[2];
sx q[2];
rz(1.549364) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0778724) q[1];
sx q[1];
rz(-1.93396) q[1];
sx q[1];
rz(2.5329464) q[1];
rz(-pi) q[2];
rz(1.767166) q[3];
sx q[3];
rz(-1.0305627) q[3];
sx q[3];
rz(-1.9062717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3123902) q[2];
sx q[2];
rz(-0.94909334) q[2];
sx q[2];
rz(0.20425805) q[2];
rz(-1.9355109) q[3];
sx q[3];
rz(-1.5217425) q[3];
sx q[3];
rz(2.9061785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4181353) q[0];
sx q[0];
rz(-1.3292987) q[0];
sx q[0];
rz(1.4468505) q[0];
rz(1.8824668) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(2.4553305) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0745569) q[0];
sx q[0];
rz(-2.0445163) q[0];
sx q[0];
rz(-2.0635598) q[0];
rz(-0.83001901) q[2];
sx q[2];
rz(-1.1168715) q[2];
sx q[2];
rz(-1.9759076) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.9069179) q[1];
sx q[1];
rz(-1.385681) q[1];
sx q[1];
rz(-0.075637416) q[1];
rz(-pi) q[2];
rz(-1.8920184) q[3];
sx q[3];
rz(-1.1676844) q[3];
sx q[3];
rz(1.5954799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.69616047) q[2];
sx q[2];
rz(-1.7636718) q[2];
sx q[2];
rz(-0.0017722842) q[2];
rz(-0.56162515) q[3];
sx q[3];
rz(-2.2300945) q[3];
sx q[3];
rz(1.6368438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5381662) q[0];
sx q[0];
rz(-0.68646938) q[0];
sx q[0];
rz(-1.6954533) q[0];
rz(0.7810477) q[1];
sx q[1];
rz(-1.8361517) q[1];
sx q[1];
rz(-1.5015645) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9672464) q[0];
sx q[0];
rz(-1.5513199) q[0];
sx q[0];
rz(2.0635701) q[0];
x q[1];
rz(-1.5090452) q[2];
sx q[2];
rz(-1.5442863) q[2];
sx q[2];
rz(0.74552958) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6540263) q[1];
sx q[1];
rz(-1.2441725) q[1];
sx q[1];
rz(-1.6649654) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99174188) q[3];
sx q[3];
rz(-1.0240882) q[3];
sx q[3];
rz(-2.9255097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.35187307) q[2];
sx q[2];
rz(-1.3616273) q[2];
sx q[2];
rz(-1.3191351) q[2];
rz(1.9296648) q[3];
sx q[3];
rz(-1.8550248) q[3];
sx q[3];
rz(0.31931988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33655745) q[0];
sx q[0];
rz(-0.55258495) q[0];
sx q[0];
rz(1.9375027) q[0];
rz(2.7583292) q[1];
sx q[1];
rz(-0.52572322) q[1];
sx q[1];
rz(-0.35167545) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6707014) q[0];
sx q[0];
rz(-2.0235217) q[0];
sx q[0];
rz(-2.7265413) q[0];
x q[1];
rz(-0.69182379) q[2];
sx q[2];
rz(-1.3335388) q[2];
sx q[2];
rz(-1.1292063) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.069236) q[1];
sx q[1];
rz(-1.623739) q[1];
sx q[1];
rz(-2.3318021) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5185235) q[3];
sx q[3];
rz(-1.6802603) q[3];
sx q[3];
rz(1.048552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.3433156) q[2];
sx q[2];
rz(-1.1077935) q[2];
sx q[2];
rz(-1.2822255) q[2];
rz(-1.4964237) q[3];
sx q[3];
rz(-1.5346425) q[3];
sx q[3];
rz(-2.0857247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4984109) q[0];
sx q[0];
rz(-1.8739941) q[0];
sx q[0];
rz(2.9472651) q[0];
rz(1.0378029) q[1];
sx q[1];
rz(-0.56832814) q[1];
sx q[1];
rz(-1.0338354) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68034222) q[0];
sx q[0];
rz(-0.67078062) q[0];
sx q[0];
rz(1.3602507) q[0];
rz(-pi) q[1];
rz(-2.6331484) q[2];
sx q[2];
rz(-2.2061081) q[2];
sx q[2];
rz(0.19257643) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6453213) q[1];
sx q[1];
rz(-1.4693345) q[1];
sx q[1];
rz(-0.99781499) q[1];
rz(2.4284027) q[3];
sx q[3];
rz(-1.4048368) q[3];
sx q[3];
rz(-0.99456577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0795435) q[2];
sx q[2];
rz(-0.94576183) q[2];
sx q[2];
rz(-2.5058084) q[2];
rz(0.27030269) q[3];
sx q[3];
rz(-2.342194) q[3];
sx q[3];
rz(1.5283782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(1.6939659) q[0];
sx q[0];
rz(-1.3128558) q[0];
sx q[0];
rz(-2.0679612) q[0];
rz(1.4355961) q[1];
sx q[1];
rz(-1.5789079) q[1];
sx q[1];
rz(0.78067738) q[1];
rz(2.1736017) q[2];
sx q[2];
rz(-1.597076) q[2];
sx q[2];
rz(1.1521641) q[2];
rz(1.2407606) q[3];
sx q[3];
rz(-1.6374554) q[3];
sx q[3];
rz(-1.0367254) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];