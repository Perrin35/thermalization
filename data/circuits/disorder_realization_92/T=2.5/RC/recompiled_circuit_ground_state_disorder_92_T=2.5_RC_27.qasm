OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.8808682) q[0];
sx q[0];
rz(3.4178419) q[0];
sx q[0];
rz(8.5810241) q[0];
rz(1.002797) q[1];
sx q[1];
rz(-2.306814) q[1];
sx q[1];
rz(0.53258449) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40520378) q[0];
sx q[0];
rz(-1.7199992) q[0];
sx q[0];
rz(-0.91513855) q[0];
x q[1];
rz(-1.4307099) q[2];
sx q[2];
rz(-1.0647961) q[2];
sx q[2];
rz(1.2337453) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0805894) q[1];
sx q[1];
rz(-2.4312651) q[1];
sx q[1];
rz(-1.3377473) q[1];
rz(-pi) q[2];
rz(-2.769773) q[3];
sx q[3];
rz(-1.9247247) q[3];
sx q[3];
rz(-1.33385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3089932) q[2];
sx q[2];
rz(-1.4048046) q[2];
sx q[2];
rz(-0.12283202) q[2];
rz(3.099856) q[3];
sx q[3];
rz(-1.5580956) q[3];
sx q[3];
rz(-2.6532555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9915344) q[0];
sx q[0];
rz(-2.8128862) q[0];
sx q[0];
rz(-2.026189) q[0];
rz(-0.54488048) q[1];
sx q[1];
rz(-1.7439525) q[1];
sx q[1];
rz(-0.56745183) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7780964) q[0];
sx q[0];
rz(-1.5370742) q[0];
sx q[0];
rz(3.0959227) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5390776) q[2];
sx q[2];
rz(-2.3467772) q[2];
sx q[2];
rz(2.1057801) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2360252) q[1];
sx q[1];
rz(-1.6710588) q[1];
sx q[1];
rz(0.57159337) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53506923) q[3];
sx q[3];
rz(-2.5949083) q[3];
sx q[3];
rz(-3.1233623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6313717) q[2];
sx q[2];
rz(-0.12237445) q[2];
sx q[2];
rz(-3.0369634) q[2];
rz(2.6369324) q[3];
sx q[3];
rz(-1.3480836) q[3];
sx q[3];
rz(0.73205718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-2.0013244) q[0];
sx q[0];
rz(-2.9502385) q[0];
sx q[0];
rz(1.6562756) q[0];
rz(-2.5775919) q[1];
sx q[1];
rz(-1.0537078) q[1];
sx q[1];
rz(-0.82411134) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018560713) q[0];
sx q[0];
rz(-2.0981952) q[0];
sx q[0];
rz(0.19077905) q[0];
rz(-2.3468909) q[2];
sx q[2];
rz(-1.1247171) q[2];
sx q[2];
rz(2.5925555) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.50269404) q[1];
sx q[1];
rz(-0.47595176) q[1];
sx q[1];
rz(-2.6902728) q[1];
rz(-pi) q[2];
rz(1.885627) q[3];
sx q[3];
rz(-0.22284976) q[3];
sx q[3];
rz(0.94161445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.092502681) q[2];
sx q[2];
rz(-0.29128942) q[2];
sx q[2];
rz(-1.5352486) q[2];
rz(2.2384079) q[3];
sx q[3];
rz(-1.817037) q[3];
sx q[3];
rz(-0.42455348) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9400738) q[0];
sx q[0];
rz(-0.98648447) q[0];
sx q[0];
rz(2.7281813) q[0];
rz(0.19551936) q[1];
sx q[1];
rz(-2.1045411) q[1];
sx q[1];
rz(1.6874541) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66178807) q[0];
sx q[0];
rz(-0.88484064) q[0];
sx q[0];
rz(0.16459008) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7614582) q[2];
sx q[2];
rz(-1.2209998) q[2];
sx q[2];
rz(1.7226532) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.11024347) q[1];
sx q[1];
rz(-1.8520466) q[1];
sx q[1];
rz(1.0737277) q[1];
rz(-pi) q[2];
rz(1.7735673) q[3];
sx q[3];
rz(-0.88651171) q[3];
sx q[3];
rz(0.87316421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2426408) q[2];
sx q[2];
rz(-2.5256038) q[2];
sx q[2];
rz(-1.2658489) q[2];
rz(-2.2856581) q[3];
sx q[3];
rz(-1.3866235) q[3];
sx q[3];
rz(-1.2148733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-0.6619381) q[0];
sx q[0];
rz(-0.18331535) q[0];
sx q[0];
rz(-1.7228175) q[0];
rz(-1.4310369) q[1];
sx q[1];
rz(-1.6173671) q[1];
sx q[1];
rz(0.72351825) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54004266) q[0];
sx q[0];
rz(-1.3941843) q[0];
sx q[0];
rz(-2.5536182) q[0];
rz(-2.3695643) q[2];
sx q[2];
rz(-2.2197422) q[2];
sx q[2];
rz(1.8963501) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11746487) q[1];
sx q[1];
rz(-0.70677838) q[1];
sx q[1];
rz(-1.5670845) q[1];
rz(0.35227065) q[3];
sx q[3];
rz(-1.5639502) q[3];
sx q[3];
rz(-0.37876836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4414703) q[2];
sx q[2];
rz(-1.0157601) q[2];
sx q[2];
rz(0.36386841) q[2];
rz(-1.8338592) q[3];
sx q[3];
rz(-1.6677083) q[3];
sx q[3];
rz(1.3522805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.42809197) q[0];
sx q[0];
rz(-2.2245753) q[0];
sx q[0];
rz(-2.2416903) q[0];
rz(-1.046754) q[1];
sx q[1];
rz(-2.7622107) q[1];
sx q[1];
rz(-2.5894763) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4259767) q[0];
sx q[0];
rz(-1.5425342) q[0];
sx q[0];
rz(0.33833276) q[0];
x q[1];
rz(-0.70163597) q[2];
sx q[2];
rz(-2.0011099) q[2];
sx q[2];
rz(2.9308191) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0585409) q[1];
sx q[1];
rz(-1.4803491) q[1];
sx q[1];
rz(1.3988711) q[1];
x q[2];
rz(2.5512929) q[3];
sx q[3];
rz(-1.5029782) q[3];
sx q[3];
rz(-0.53896157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3209352) q[2];
sx q[2];
rz(-0.67639095) q[2];
sx q[2];
rz(-0.19860849) q[2];
rz(1.1292388) q[3];
sx q[3];
rz(-1.4720474) q[3];
sx q[3];
rz(-0.78013295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3405014) q[0];
sx q[0];
rz(-1.3022364) q[0];
sx q[0];
rz(2.4275725) q[0];
rz(-1.2227614) q[1];
sx q[1];
rz(-2.265265) q[1];
sx q[1];
rz(-0.27539918) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.794316) q[0];
sx q[0];
rz(-1.1755318) q[0];
sx q[0];
rz(-2.9121132) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95850079) q[2];
sx q[2];
rz(-1.3103518) q[2];
sx q[2];
rz(0.97439235) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.279567) q[1];
sx q[1];
rz(-1.0890402) q[1];
sx q[1];
rz(-0.53558366) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.015206788) q[3];
sx q[3];
rz(-0.72532024) q[3];
sx q[3];
rz(-2.050658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3068646) q[2];
sx q[2];
rz(-1.6708259) q[2];
sx q[2];
rz(1.6738221) q[2];
rz(1.62014) q[3];
sx q[3];
rz(-0.68632564) q[3];
sx q[3];
rz(0.93792382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15243212) q[0];
sx q[0];
rz(-0.91355649) q[0];
sx q[0];
rz(0.96617019) q[0];
rz(0.78978157) q[1];
sx q[1];
rz(-0.25895324) q[1];
sx q[1];
rz(-2.1678179) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1137266) q[0];
sx q[0];
rz(-1.6902414) q[0];
sx q[0];
rz(1.8177423) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70232161) q[2];
sx q[2];
rz(-0.53793814) q[2];
sx q[2];
rz(0.59662102) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3605347) q[1];
sx q[1];
rz(-1.2670749) q[1];
sx q[1];
rz(-2.2045434) q[1];
rz(-3.0816541) q[3];
sx q[3];
rz(-1.0956683) q[3];
sx q[3];
rz(-2.5157473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3070273) q[2];
sx q[2];
rz(-1.7119188) q[2];
sx q[2];
rz(-0.60565051) q[2];
rz(2.0511131) q[3];
sx q[3];
rz(-2.8993789) q[3];
sx q[3];
rz(-0.098701326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7379446) q[0];
sx q[0];
rz(-0.048796766) q[0];
sx q[0];
rz(-0.57156372) q[0];
rz(2.7538708) q[1];
sx q[1];
rz(-2.3523836) q[1];
sx q[1];
rz(0.8943843) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26167241) q[0];
sx q[0];
rz(-2.399647) q[0];
sx q[0];
rz(0.80961734) q[0];
x q[1];
rz(1.0036976) q[2];
sx q[2];
rz(-1.4722479) q[2];
sx q[2];
rz(-0.71426094) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.91796658) q[1];
sx q[1];
rz(-1.8119748) q[1];
sx q[1];
rz(1.1583223) q[1];
rz(1.649141) q[3];
sx q[3];
rz(-2.0781029) q[3];
sx q[3];
rz(2.0149734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1423219) q[2];
sx q[2];
rz(-2.2597376) q[2];
sx q[2];
rz(-0.87232653) q[2];
rz(-0.074660389) q[3];
sx q[3];
rz(-1.9118237) q[3];
sx q[3];
rz(-2.5653896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1416624) q[0];
sx q[0];
rz(-0.44121989) q[0];
sx q[0];
rz(-0.73750752) q[0];
rz(-2.4420338) q[1];
sx q[1];
rz(-2.12205) q[1];
sx q[1];
rz(-2.2950744) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93447157) q[0];
sx q[0];
rz(-1.1729329) q[0];
sx q[0];
rz(-2.4475696) q[0];
rz(1.9969127) q[2];
sx q[2];
rz(-1.6175381) q[2];
sx q[2];
rz(-0.5668723) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.10900765) q[1];
sx q[1];
rz(-1.4334588) q[1];
sx q[1];
rz(2.0020383) q[1];
rz(-0.41937866) q[3];
sx q[3];
rz(-0.43361317) q[3];
sx q[3];
rz(2.8948642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.49125853) q[2];
sx q[2];
rz(-1.9776055) q[2];
sx q[2];
rz(-0.35438806) q[2];
rz(0.77704159) q[3];
sx q[3];
rz(-0.6207501) q[3];
sx q[3];
rz(1.0907772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(-2.1914094) q[0];
sx q[0];
rz(-1.2744899) q[0];
sx q[0];
rz(0.50344678) q[0];
rz(1.3632111) q[1];
sx q[1];
rz(-2.6111205) q[1];
sx q[1];
rz(3.0725239) q[1];
rz(-3.0341078) q[2];
sx q[2];
rz(-0.86089118) q[2];
sx q[2];
rz(-0.64313342) q[2];
rz(0.94181413) q[3];
sx q[3];
rz(-2.5968646) q[3];
sx q[3];
rz(-1.4379848) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
