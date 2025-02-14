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
rz(-0.80225575) q[0];
sx q[0];
rz(1.383961) q[0];
sx q[0];
rz(11.297754) q[0];
rz(0.4624548) q[1];
sx q[1];
rz(6.3422536) q[1];
sx q[1];
rz(8.0191945) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6574638) q[0];
sx q[0];
rz(-1.1384369) q[0];
sx q[0];
rz(2.0922489) q[0];
rz(-pi) q[1];
rz(-0.67063801) q[2];
sx q[2];
rz(-1.7522289) q[2];
sx q[2];
rz(-1.8701613) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.51912145) q[1];
sx q[1];
rz(-2.2833385) q[1];
sx q[1];
rz(2.126241) q[1];
rz(2.9146092) q[3];
sx q[3];
rz(-0.90714083) q[3];
sx q[3];
rz(-3.0945145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1543701) q[2];
sx q[2];
rz(-1.1673678) q[2];
sx q[2];
rz(2.3517189) q[2];
rz(3.1176944) q[3];
sx q[3];
rz(-1.3393211) q[3];
sx q[3];
rz(2.6051615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8928878) q[0];
sx q[0];
rz(-1.6510115) q[0];
sx q[0];
rz(-1.8875246) q[0];
rz(1.4887384) q[1];
sx q[1];
rz(-0.87569991) q[1];
sx q[1];
rz(0.48315963) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41071677) q[0];
sx q[0];
rz(-2.9456101) q[0];
sx q[0];
rz(1.8394801) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0873726) q[2];
sx q[2];
rz(-0.90460515) q[2];
sx q[2];
rz(1.7214799) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.28774481) q[1];
sx q[1];
rz(-1.8327296) q[1];
sx q[1];
rz(0.37516428) q[1];
rz(1.7464569) q[3];
sx q[3];
rz(-1.2845728) q[3];
sx q[3];
rz(-1.3940879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.88863215) q[2];
sx q[2];
rz(-1.4560207) q[2];
sx q[2];
rz(1.6271094) q[2];
rz(3.0857981) q[3];
sx q[3];
rz(-1.5443708) q[3];
sx q[3];
rz(-1.4896711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(2.2273939) q[0];
sx q[0];
rz(-2.6970503) q[0];
sx q[0];
rz(3.0134873) q[0];
rz(0.57614342) q[1];
sx q[1];
rz(-2.4100401) q[1];
sx q[1];
rz(0.76549706) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.573581) q[0];
sx q[0];
rz(-0.67493248) q[0];
sx q[0];
rz(0.57008596) q[0];
x q[1];
rz(-2.4548495) q[2];
sx q[2];
rz(-1.8675065) q[2];
sx q[2];
rz(0.77881694) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.004359) q[1];
sx q[1];
rz(-0.69789125) q[1];
sx q[1];
rz(-0.76142074) q[1];
x q[2];
rz(-1.5452905) q[3];
sx q[3];
rz(-1.4273941) q[3];
sx q[3];
rz(2.0633782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3196044) q[2];
sx q[2];
rz(-0.87279785) q[2];
sx q[2];
rz(2.2785211) q[2];
rz(-0.34267628) q[3];
sx q[3];
rz(-2.7211012) q[3];
sx q[3];
rz(-1.4884523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89219013) q[0];
sx q[0];
rz(-2.1718195) q[0];
sx q[0];
rz(-0.48200193) q[0];
rz(-0.30872289) q[1];
sx q[1];
rz(-0.32544193) q[1];
sx q[1];
rz(-2.5293005) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9252602) q[0];
sx q[0];
rz(-1.060492) q[0];
sx q[0];
rz(0.55522529) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6545534) q[2];
sx q[2];
rz(-1.5847907) q[2];
sx q[2];
rz(2.9521041) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.051381342) q[1];
sx q[1];
rz(-1.6156229) q[1];
sx q[1];
rz(-0.96446891) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5641698) q[3];
sx q[3];
rz(-1.7200791) q[3];
sx q[3];
rz(-1.7501259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.65773949) q[2];
sx q[2];
rz(-2.190399) q[2];
sx q[2];
rz(1.308002) q[2];
rz(-2.7075503) q[3];
sx q[3];
rz(-1.3515892) q[3];
sx q[3];
rz(-1.2691809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16172116) q[0];
sx q[0];
rz(-1.5525818) q[0];
sx q[0];
rz(0.27444926) q[0];
rz(1.5212003) q[1];
sx q[1];
rz(-1.043964) q[1];
sx q[1];
rz(1.8843947) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8306106) q[0];
sx q[0];
rz(-2.794696) q[0];
sx q[0];
rz(0.78126379) q[0];
rz(-pi) q[1];
rz(1.0368226) q[2];
sx q[2];
rz(-1.8920915) q[2];
sx q[2];
rz(-0.23201135) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8179743) q[1];
sx q[1];
rz(-2.5516641) q[1];
sx q[1];
rz(1.9650311) q[1];
rz(-pi) q[2];
rz(-0.86602314) q[3];
sx q[3];
rz(-1.8187722) q[3];
sx q[3];
rz(0.06125227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3955128) q[2];
sx q[2];
rz(-1.6359676) q[2];
sx q[2];
rz(-2.5858322) q[2];
rz(2.3505576) q[3];
sx q[3];
rz(-1.0020703) q[3];
sx q[3];
rz(1.6382431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.898734) q[0];
sx q[0];
rz(-1.4396311) q[0];
sx q[0];
rz(-2.4295501) q[0];
rz(-2.5391319) q[1];
sx q[1];
rz(-2.7280877) q[1];
sx q[1];
rz(-3.0625524) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37385434) q[0];
sx q[0];
rz(-1.615623) q[0];
sx q[0];
rz(1.5064513) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4858133) q[2];
sx q[2];
rz(-1.4956317) q[2];
sx q[2];
rz(-2.4369881) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.5309294) q[1];
sx q[1];
rz(-0.77745118) q[1];
sx q[1];
rz(2.8270097) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32460605) q[3];
sx q[3];
rz(-0.31346009) q[3];
sx q[3];
rz(2.4544668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.041212335) q[2];
sx q[2];
rz(-2.2104287) q[2];
sx q[2];
rz(2.1616914) q[2];
rz(-1.529871) q[3];
sx q[3];
rz(-1.2414705) q[3];
sx q[3];
rz(2.8821168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74705446) q[0];
sx q[0];
rz(-1.2815481) q[0];
sx q[0];
rz(3.0406612) q[0];
rz(-2.586567) q[1];
sx q[1];
rz(-0.9340159) q[1];
sx q[1];
rz(-1.8468044) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0880343) q[0];
sx q[0];
rz(-2.4248666) q[0];
sx q[0];
rz(1.2117366) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2764658) q[2];
sx q[2];
rz(-1.571725) q[2];
sx q[2];
rz(-1.1565191) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.87161049) q[1];
sx q[1];
rz(-1.9967604) q[1];
sx q[1];
rz(1.4130089) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8517868) q[3];
sx q[3];
rz(-1.13124) q[3];
sx q[3];
rz(1.2729147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1780221) q[2];
sx q[2];
rz(-0.87248412) q[2];
sx q[2];
rz(-2.8998609) q[2];
rz(2.669615) q[3];
sx q[3];
rz(-0.21502544) q[3];
sx q[3];
rz(-1.5836466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(0.80411512) q[0];
sx q[0];
rz(-0.98144704) q[0];
sx q[0];
rz(2.5313983) q[0];
rz(0.21774165) q[1];
sx q[1];
rz(-2.1861031) q[1];
sx q[1];
rz(-2.311923) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1581081) q[0];
sx q[0];
rz(-0.89644855) q[0];
sx q[0];
rz(-0.83263065) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2118633) q[2];
sx q[2];
rz(-1.8808489) q[2];
sx q[2];
rz(-0.25676766) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9431994) q[1];
sx q[1];
rz(-1.3855318) q[1];
sx q[1];
rz(-2.7279502) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85932087) q[3];
sx q[3];
rz(-1.4830657) q[3];
sx q[3];
rz(0.52857196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3024451) q[2];
sx q[2];
rz(-1.382901) q[2];
sx q[2];
rz(-1.1350606) q[2];
rz(-0.48842397) q[3];
sx q[3];
rz(-1.6596183) q[3];
sx q[3];
rz(1.9058913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9957073) q[0];
sx q[0];
rz(-1.1111525) q[0];
sx q[0];
rz(2.1671894) q[0];
rz(-0.79565945) q[1];
sx q[1];
rz(-0.37745044) q[1];
sx q[1];
rz(-2.5239351) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97516135) q[0];
sx q[0];
rz(-1.5733061) q[0];
sx q[0];
rz(3.0581492) q[0];
rz(-1.9699924) q[2];
sx q[2];
rz(-2.0621057) q[2];
sx q[2];
rz(-2.0727061) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1682081) q[1];
sx q[1];
rz(-2.1323556) q[1];
sx q[1];
rz(1.8902768) q[1];
x q[2];
rz(-1.0997201) q[3];
sx q[3];
rz(-1.866939) q[3];
sx q[3];
rz(-2.5961034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1254897) q[2];
sx q[2];
rz(-1.2697271) q[2];
sx q[2];
rz(-1.7529091) q[2];
rz(1.8359418) q[3];
sx q[3];
rz(-0.7235705) q[3];
sx q[3];
rz(1.2828264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1284803) q[0];
sx q[0];
rz(-2.4025669) q[0];
sx q[0];
rz(-0.59984961) q[0];
rz(-2.5685617) q[1];
sx q[1];
rz(-2.532798) q[1];
sx q[1];
rz(-1.0240239) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065627873) q[0];
sx q[0];
rz(-0.5813404) q[0];
sx q[0];
rz(2.3443048) q[0];
x q[1];
rz(-3.0863831) q[2];
sx q[2];
rz(-2.488236) q[2];
sx q[2];
rz(-0.094815985) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.21067317) q[1];
sx q[1];
rz(-1.966488) q[1];
sx q[1];
rz(-1.4798505) q[1];
rz(-1.4198331) q[3];
sx q[3];
rz(-1.9486041) q[3];
sx q[3];
rz(-1.0491766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2188501) q[2];
sx q[2];
rz(-1.7538193) q[2];
sx q[2];
rz(-0.77524033) q[2];
rz(1.5357337) q[3];
sx q[3];
rz(-1.1865059) q[3];
sx q[3];
rz(1.8839914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53104644) q[0];
sx q[0];
rz(-2.1903867) q[0];
sx q[0];
rz(1.2280986) q[0];
rz(1.3950521) q[1];
sx q[1];
rz(-1.6680622) q[1];
sx q[1];
rz(-0.39573085) q[1];
rz(-0.45356815) q[2];
sx q[2];
rz(-1.8487052) q[2];
sx q[2];
rz(-2.1089274) q[2];
rz(3.1093194) q[3];
sx q[3];
rz(-2.3393031) q[3];
sx q[3];
rz(2.5325251) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
