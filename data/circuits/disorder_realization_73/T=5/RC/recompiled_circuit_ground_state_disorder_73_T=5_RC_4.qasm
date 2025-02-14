OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0176528) q[0];
sx q[0];
rz(-2.2247563) q[0];
sx q[0];
rz(0.43489781) q[0];
rz(1.5974367) q[1];
sx q[1];
rz(3.6606001) q[1];
sx q[1];
rz(11.984348) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30729025) q[0];
sx q[0];
rz(-2.8367762) q[0];
sx q[0];
rz(2.1075664) q[0];
x q[1];
rz(-2.2088654) q[2];
sx q[2];
rz(-1.1865532) q[2];
sx q[2];
rz(-2.8926433) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.21596244) q[1];
sx q[1];
rz(-1.1396038) q[1];
sx q[1];
rz(2.870382) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9635756) q[3];
sx q[3];
rz(-1.2937163) q[3];
sx q[3];
rz(1.8711562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.55328289) q[2];
sx q[2];
rz(-1.8841691) q[2];
sx q[2];
rz(2.8498939) q[2];
rz(2.6907673) q[3];
sx q[3];
rz(-2.8901633) q[3];
sx q[3];
rz(-2.5341471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70158231) q[0];
sx q[0];
rz(-2.1511183) q[0];
sx q[0];
rz(2.0462346) q[0];
rz(-1.9502684) q[1];
sx q[1];
rz(-2.2215863) q[1];
sx q[1];
rz(0.90257588) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1987225) q[0];
sx q[0];
rz(-0.96317055) q[0];
sx q[0];
rz(-2.3234419) q[0];
rz(0.90229374) q[2];
sx q[2];
rz(-1.3375139) q[2];
sx q[2];
rz(-1.4177314) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2486156) q[1];
sx q[1];
rz(-1.4505523) q[1];
sx q[1];
rz(-0.29386947) q[1];
rz(0.9489305) q[3];
sx q[3];
rz(-1.9919259) q[3];
sx q[3];
rz(1.3768743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.1878745) q[2];
sx q[2];
rz(-2.1672921) q[2];
sx q[2];
rz(-0.11889674) q[2];
rz(2.4925354) q[3];
sx q[3];
rz(-0.18638149) q[3];
sx q[3];
rz(-1.4547179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.651799) q[0];
sx q[0];
rz(-1.8875903) q[0];
sx q[0];
rz(-0.5157665) q[0];
rz(2.0491397) q[1];
sx q[1];
rz(-0.40320435) q[1];
sx q[1];
rz(-1.2752424) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9495205) q[0];
sx q[0];
rz(-1.6123591) q[0];
sx q[0];
rz(-1.5526718) q[0];
rz(0.31351201) q[2];
sx q[2];
rz(-1.3517153) q[2];
sx q[2];
rz(-0.30922019) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0787492) q[1];
sx q[1];
rz(-0.8991407) q[1];
sx q[1];
rz(-0.25154227) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10026513) q[3];
sx q[3];
rz(-2.3449922) q[3];
sx q[3];
rz(-1.0880566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0755997) q[2];
sx q[2];
rz(-1.8241355) q[2];
sx q[2];
rz(-1.9817748) q[2];
rz(-2.6521111) q[3];
sx q[3];
rz(-1.5882086) q[3];
sx q[3];
rz(-0.71973962) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2751665) q[0];
sx q[0];
rz(-2.8186099) q[0];
sx q[0];
rz(2.6599356) q[0];
rz(-0.9306759) q[1];
sx q[1];
rz(-1.8447256) q[1];
sx q[1];
rz(-0.01893386) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7304283) q[0];
sx q[0];
rz(-1.5207964) q[0];
sx q[0];
rz(3.1206112) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2332814) q[2];
sx q[2];
rz(-1.5358972) q[2];
sx q[2];
rz(2.1878583) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1065797) q[1];
sx q[1];
rz(-2.0845319) q[1];
sx q[1];
rz(-2.4519743) q[1];
rz(2.3871759) q[3];
sx q[3];
rz(-2.2212611) q[3];
sx q[3];
rz(-2.2353896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9229752) q[2];
sx q[2];
rz(-2.7851892) q[2];
sx q[2];
rz(2.3962928) q[2];
rz(-2.7636012) q[3];
sx q[3];
rz(-1.1443006) q[3];
sx q[3];
rz(2.2530344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22076631) q[0];
sx q[0];
rz(-2.8296318) q[0];
sx q[0];
rz(1.9507324) q[0];
rz(-0.4885172) q[1];
sx q[1];
rz(-2.2256475) q[1];
sx q[1];
rz(-0.8078422) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4427748) q[0];
sx q[0];
rz(-1.6614058) q[0];
sx q[0];
rz(-1.242799) q[0];
x q[1];
rz(-1.135181) q[2];
sx q[2];
rz(-0.50335889) q[2];
sx q[2];
rz(1.8653009) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78524974) q[1];
sx q[1];
rz(-1.4828353) q[1];
sx q[1];
rz(-0.45694074) q[1];
x q[2];
rz(-0.2000974) q[3];
sx q[3];
rz(-0.95733023) q[3];
sx q[3];
rz(0.6544978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0025582) q[2];
sx q[2];
rz(-2.148874) q[2];
sx q[2];
rz(-2.9774418) q[2];
rz(2.3955591) q[3];
sx q[3];
rz(-1.7697216) q[3];
sx q[3];
rz(0.073808864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.98110759) q[0];
sx q[0];
rz(-1.2657413) q[0];
sx q[0];
rz(-2.4142081) q[0];
rz(-2.6630867) q[1];
sx q[1];
rz(-1.6886657) q[1];
sx q[1];
rz(0.7368288) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3664624) q[0];
sx q[0];
rz(-1.4331761) q[0];
sx q[0];
rz(-3.0855975) q[0];
x q[1];
rz(1.1428035) q[2];
sx q[2];
rz(-0.30747947) q[2];
sx q[2];
rz(2.3124419) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3697065) q[1];
sx q[1];
rz(-0.73644887) q[1];
sx q[1];
rz(2.1320029) q[1];
rz(-1.9584452) q[3];
sx q[3];
rz(-1.5836925) q[3];
sx q[3];
rz(1.8356334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.14787802) q[2];
sx q[2];
rz(-2.1174049) q[2];
sx q[2];
rz(2.4564339) q[2];
rz(2.1327175) q[3];
sx q[3];
rz(-2.4515371) q[3];
sx q[3];
rz(2.8907997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.93609) q[0];
sx q[0];
rz(-0.093955366) q[0];
sx q[0];
rz(2.0482546) q[0];
rz(-0.023748485) q[1];
sx q[1];
rz(-0.7370342) q[1];
sx q[1];
rz(2.645983) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0422875) q[0];
sx q[0];
rz(-1.5047502) q[0];
sx q[0];
rz(1.5586067) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43735403) q[2];
sx q[2];
rz(-1.3833356) q[2];
sx q[2];
rz(-2.1700493) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.88265282) q[1];
sx q[1];
rz(-1.7523089) q[1];
sx q[1];
rz(2.7869999) q[1];
rz(-2.7109259) q[3];
sx q[3];
rz(-2.0206631) q[3];
sx q[3];
rz(-1.7783742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2451943) q[2];
sx q[2];
rz(-2.3485025) q[2];
sx q[2];
rz(-2.8450361) q[2];
rz(-0.5101997) q[3];
sx q[3];
rz(-1.5254131) q[3];
sx q[3];
rz(-0.27142522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1475913) q[0];
sx q[0];
rz(-2.1864102) q[0];
sx q[0];
rz(-1.9872794) q[0];
rz(1.1555903) q[1];
sx q[1];
rz(-0.65735936) q[1];
sx q[1];
rz(-1.4591699) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71135581) q[0];
sx q[0];
rz(-2.2855496) q[0];
sx q[0];
rz(-1.6011333) q[0];
rz(-3.1404678) q[2];
sx q[2];
rz(-0.84648593) q[2];
sx q[2];
rz(1.6206737) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3633948) q[1];
sx q[1];
rz(-2.6043335) q[1];
sx q[1];
rz(-1.926653) q[1];
rz(-0.67787637) q[3];
sx q[3];
rz(-2.3837485) q[3];
sx q[3];
rz(2.0539396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.83361721) q[2];
sx q[2];
rz(-0.54215446) q[2];
sx q[2];
rz(-1.7657492) q[2];
rz(-0.086325072) q[3];
sx q[3];
rz(-0.83266801) q[3];
sx q[3];
rz(-1.8858006) q[3];
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
rz(1.9432705) q[0];
sx q[0];
rz(-0.69381303) q[0];
sx q[0];
rz(2.2721403) q[0];
rz(-0.72932875) q[1];
sx q[1];
rz(-0.84356934) q[1];
sx q[1];
rz(-2.9076911) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9651523) q[0];
sx q[0];
rz(-1.5230706) q[0];
sx q[0];
rz(-0.41712572) q[0];
rz(0.21017615) q[2];
sx q[2];
rz(-1.4769722) q[2];
sx q[2];
rz(1.891234) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0879957) q[1];
sx q[1];
rz(-1.1302179) q[1];
sx q[1];
rz(1.7338899) q[1];
rz(-pi) q[2];
rz(-1.9332658) q[3];
sx q[3];
rz(-2.3911282) q[3];
sx q[3];
rz(-3.124231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.087661155) q[2];
sx q[2];
rz(-1.4940741) q[2];
sx q[2];
rz(1.1762478) q[2];
rz(0.67115274) q[3];
sx q[3];
rz(-1.8920369) q[3];
sx q[3];
rz(-1.6254856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7072356) q[0];
sx q[0];
rz(-2.2800627) q[0];
sx q[0];
rz(2.0980515) q[0];
rz(3.0442944) q[1];
sx q[1];
rz(-2.0847335) q[1];
sx q[1];
rz(1.1204488) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53172511) q[0];
sx q[0];
rz(-1.5104482) q[0];
sx q[0];
rz(-2.1670695) q[0];
rz(-2.3612622) q[2];
sx q[2];
rz(-2.188211) q[2];
sx q[2];
rz(-1.8898026) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6594362) q[1];
sx q[1];
rz(-1.5220571) q[1];
sx q[1];
rz(2.963889) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2213732) q[3];
sx q[3];
rz(-0.50822483) q[3];
sx q[3];
rz(1.0168545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6843159) q[2];
sx q[2];
rz(-2.5280759) q[2];
sx q[2];
rz(1.3555917) q[2];
rz(-1.5654303) q[3];
sx q[3];
rz(-2.0578945) q[3];
sx q[3];
rz(-2.1255597) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9044357) q[0];
sx q[0];
rz(-2.5813527) q[0];
sx q[0];
rz(-0.78975633) q[0];
rz(0.65658983) q[1];
sx q[1];
rz(-1.1142535) q[1];
sx q[1];
rz(2.5744892) q[1];
rz(2.1661027) q[2];
sx q[2];
rz(-2.3885623) q[2];
sx q[2];
rz(1.9293712) q[2];
rz(-1.7651032) q[3];
sx q[3];
rz(-1.2855218) q[3];
sx q[3];
rz(2.2900477) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
