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
rz(-2.4177457) q[0];
sx q[0];
rz(-1.9404193) q[0];
sx q[0];
rz(2.6475651) q[0];
rz(1.5552893) q[1];
sx q[1];
rz(2.2145693) q[1];
sx q[1];
rz(8.5742843) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5587529) q[0];
sx q[0];
rz(-2.2338998) q[0];
sx q[0];
rz(-0.21132529) q[0];
rz(-0.57588864) q[2];
sx q[2];
rz(-0.73969361) q[2];
sx q[2];
rz(2.0832555) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.96877484) q[1];
sx q[1];
rz(-2.7451395) q[1];
sx q[1];
rz(0.60703599) q[1];
x q[2];
rz(-0.36601074) q[3];
sx q[3];
rz(-1.5029534) q[3];
sx q[3];
rz(-2.3930166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4832619) q[2];
sx q[2];
rz(-0.61507812) q[2];
sx q[2];
rz(2.1966546) q[2];
rz(-3.1350709) q[3];
sx q[3];
rz(-2.3801453) q[3];
sx q[3];
rz(0.25750461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1061123) q[0];
sx q[0];
rz(-2.1529614) q[0];
sx q[0];
rz(-0.60580564) q[0];
rz(0.93217355) q[1];
sx q[1];
rz(-1.4235539) q[1];
sx q[1];
rz(-0.1056284) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6208134) q[0];
sx q[0];
rz(-2.438669) q[0];
sx q[0];
rz(-0.23971324) q[0];
rz(-pi) q[1];
rz(-2.7850161) q[2];
sx q[2];
rz(-2.0322213) q[2];
sx q[2];
rz(0.8666641) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.88813284) q[1];
sx q[1];
rz(-1.8646629) q[1];
sx q[1];
rz(-0.47305686) q[1];
rz(-pi) q[2];
x q[2];
rz(0.83189957) q[3];
sx q[3];
rz(-0.53153342) q[3];
sx q[3];
rz(-0.040415045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2848844) q[2];
sx q[2];
rz(-1.363089) q[2];
sx q[2];
rz(-1.523783) q[2];
rz(-2.3273322) q[3];
sx q[3];
rz(-1.3596478) q[3];
sx q[3];
rz(-2.2566569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.043269) q[0];
sx q[0];
rz(-2.3461778) q[0];
sx q[0];
rz(1.5765618) q[0];
rz(0.99524975) q[1];
sx q[1];
rz(-2.1627656) q[1];
sx q[1];
rz(-2.3547122) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7835307) q[0];
sx q[0];
rz(-1.3710877) q[0];
sx q[0];
rz(3.003503) q[0];
rz(-pi) q[1];
rz(1.4866531) q[2];
sx q[2];
rz(-1.9168789) q[2];
sx q[2];
rz(3.0402407) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.39549024) q[1];
sx q[1];
rz(-1.8999892) q[1];
sx q[1];
rz(-1.2389061) q[1];
rz(-pi) q[2];
rz(-0.78600471) q[3];
sx q[3];
rz(-0.084298221) q[3];
sx q[3];
rz(0.5072197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3028822) q[2];
sx q[2];
rz(-1.6532712) q[2];
sx q[2];
rz(-0.6380471) q[2];
rz(2.6436515) q[3];
sx q[3];
rz(-2.0697856) q[3];
sx q[3];
rz(2.0518484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2447253) q[0];
sx q[0];
rz(-1.7843972) q[0];
sx q[0];
rz(3.0778399) q[0];
rz(-1.6925192) q[1];
sx q[1];
rz(-1.8359102) q[1];
sx q[1];
rz(1.8316899) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5638014) q[0];
sx q[0];
rz(-0.088106958) q[0];
sx q[0];
rz(-2.7035575) q[0];
x q[1];
rz(1.3856349) q[2];
sx q[2];
rz(-1.7660391) q[2];
sx q[2];
rz(2.756292) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.79015649) q[1];
sx q[1];
rz(-0.75439851) q[1];
sx q[1];
rz(0.18138563) q[1];
x q[2];
rz(-1.1442361) q[3];
sx q[3];
rz(-2.4080695) q[3];
sx q[3];
rz(-0.13338062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0706851) q[2];
sx q[2];
rz(-2.0347774) q[2];
sx q[2];
rz(3.0359388) q[2];
rz(1.9581155) q[3];
sx q[3];
rz(-0.89619291) q[3];
sx q[3];
rz(-2.4699874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79528177) q[0];
sx q[0];
rz(-2.2310937) q[0];
sx q[0];
rz(-1.7972535) q[0];
rz(-2.4686939) q[1];
sx q[1];
rz(-1.3105323) q[1];
sx q[1];
rz(-0.013462822) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1707662) q[0];
sx q[0];
rz(-1.544569) q[0];
sx q[0];
rz(-2.2607525) q[0];
rz(-pi) q[1];
rz(-3.1273753) q[2];
sx q[2];
rz(-1.1566312) q[2];
sx q[2];
rz(-2.3384475) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8103761) q[1];
sx q[1];
rz(-1.5575641) q[1];
sx q[1];
rz(2.7262336) q[1];
rz(2.2736566) q[3];
sx q[3];
rz(-1.867467) q[3];
sx q[3];
rz(2.5698667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.54231918) q[2];
sx q[2];
rz(-0.66212526) q[2];
sx q[2];
rz(-1.8947961) q[2];
rz(-3.0689734) q[3];
sx q[3];
rz(-2.9473372) q[3];
sx q[3];
rz(-0.3092002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70960629) q[0];
sx q[0];
rz(-1.1259587) q[0];
sx q[0];
rz(3.0288938) q[0];
rz(-0.20507774) q[1];
sx q[1];
rz(-1.3321184) q[1];
sx q[1];
rz(2.233706) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8255709) q[0];
sx q[0];
rz(-2.1793723) q[0];
sx q[0];
rz(0.84644239) q[0];
rz(-2.5656469) q[2];
sx q[2];
rz(-2.4436946) q[2];
sx q[2];
rz(-0.18651785) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9428633) q[1];
sx q[1];
rz(-2.2009549) q[1];
sx q[1];
rz(-1.4575926) q[1];
rz(1.7849633) q[3];
sx q[3];
rz(-2.2244144) q[3];
sx q[3];
rz(0.71398338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9997361) q[2];
sx q[2];
rz(-0.33679589) q[2];
sx q[2];
rz(3.0736308) q[2];
rz(1.9939907) q[3];
sx q[3];
rz(-1.8736898) q[3];
sx q[3];
rz(1.2609153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9481908) q[0];
sx q[0];
rz(-1.3113439) q[0];
sx q[0];
rz(-2.6780658) q[0];
rz(-2.9947128) q[1];
sx q[1];
rz(-1.2356707) q[1];
sx q[1];
rz(-0.99259496) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3459839) q[0];
sx q[0];
rz(-0.75712403) q[0];
sx q[0];
rz(0.36619314) q[0];
rz(-pi) q[1];
rz(1.2847177) q[2];
sx q[2];
rz(-1.5062638) q[2];
sx q[2];
rz(-0.71071029) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.22659644) q[1];
sx q[1];
rz(-2.2961756) q[1];
sx q[1];
rz(-0.62208773) q[1];
rz(0.73294183) q[3];
sx q[3];
rz(-1.2228726) q[3];
sx q[3];
rz(-0.80845736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.85840449) q[2];
sx q[2];
rz(-1.4294727) q[2];
sx q[2];
rz(0.44537133) q[2];
rz(-1.3238268) q[3];
sx q[3];
rz(-2.0711074) q[3];
sx q[3];
rz(-2.0557192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0367947) q[0];
sx q[0];
rz(-0.0093655149) q[0];
sx q[0];
rz(2.985756) q[0];
rz(1.8644631) q[1];
sx q[1];
rz(-1.7946449) q[1];
sx q[1];
rz(-3.1256622) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30980047) q[0];
sx q[0];
rz(-2.3224717) q[0];
sx q[0];
rz(2.8099634) q[0];
x q[1];
rz(-2.7819949) q[2];
sx q[2];
rz(-2.3162875) q[2];
sx q[2];
rz(-1.8458837) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.58305743) q[1];
sx q[1];
rz(-1.0638405) q[1];
sx q[1];
rz(0.81113775) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5032477) q[3];
sx q[3];
rz(-1.2223635) q[3];
sx q[3];
rz(-1.3909457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5126123) q[2];
sx q[2];
rz(-3.0292558) q[2];
sx q[2];
rz(-0.12574276) q[2];
rz(0.9564774) q[3];
sx q[3];
rz(-1.4230909) q[3];
sx q[3];
rz(2.8023348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9823031) q[0];
sx q[0];
rz(-0.24022261) q[0];
sx q[0];
rz(-2.6851728) q[0];
rz(1.5688815) q[1];
sx q[1];
rz(-1.0036889) q[1];
sx q[1];
rz(-2.454954) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7266885) q[0];
sx q[0];
rz(-0.81302887) q[0];
sx q[0];
rz(-0.70720478) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1988611) q[2];
sx q[2];
rz(-0.21272993) q[2];
sx q[2];
rz(-2.5036734) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9514044) q[1];
sx q[1];
rz(-1.4238289) q[1];
sx q[1];
rz(1.4495871) q[1];
rz(-0.37952642) q[3];
sx q[3];
rz(-0.18774334) q[3];
sx q[3];
rz(1.2953341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.75866428) q[2];
sx q[2];
rz(-1.1229346) q[2];
sx q[2];
rz(2.0056966) q[2];
rz(-2.1618333) q[3];
sx q[3];
rz(-1.3330678) q[3];
sx q[3];
rz(-0.69825828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6887688) q[0];
sx q[0];
rz(-1.3249506) q[0];
sx q[0];
rz(0.45595566) q[0];
rz(-3.0988354) q[1];
sx q[1];
rz(-1.9330934) q[1];
sx q[1];
rz(2.4483689) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15007818) q[0];
sx q[0];
rz(-2.743609) q[0];
sx q[0];
rz(1.8282169) q[0];
rz(-pi) q[1];
rz(-1.1666388) q[2];
sx q[2];
rz(-1.9459138) q[2];
sx q[2];
rz(-2.1585495) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7267513) q[1];
sx q[1];
rz(-1.5637239) q[1];
sx q[1];
rz(1.4727791) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59488036) q[3];
sx q[3];
rz(-0.52973807) q[3];
sx q[3];
rz(1.6486419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2263055) q[2];
sx q[2];
rz(-1.8733571) q[2];
sx q[2];
rz(2.477296) q[2];
rz(1.9281467) q[3];
sx q[3];
rz(-2.4197141) q[3];
sx q[3];
rz(2.2693995) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.450347) q[0];
sx q[0];
rz(-2.3007614) q[0];
sx q[0];
rz(1.7578516) q[0];
rz(-1.6830403) q[1];
sx q[1];
rz(-0.28266193) q[1];
sx q[1];
rz(-2.2255486) q[1];
rz(3.0677879) q[2];
sx q[2];
rz(-2.7701785) q[2];
sx q[2];
rz(-2.2899173) q[2];
rz(-2.1324329) q[3];
sx q[3];
rz(-2.1957993) q[3];
sx q[3];
rz(-1.5002863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
