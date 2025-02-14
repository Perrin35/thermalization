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
rz(0.63737386) q[0];
sx q[0];
rz(7.3399788) q[0];
sx q[0];
rz(10.196431) q[0];
rz(-0.15089384) q[1];
sx q[1];
rz(-0.3124736) q[1];
sx q[1];
rz(1.4788652) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85017289) q[0];
sx q[0];
rz(-1.2787191) q[0];
sx q[0];
rz(0.34698457) q[0];
rz(-pi) q[1];
rz(2.0059228) q[2];
sx q[2];
rz(-0.56081334) q[2];
sx q[2];
rz(-0.31080526) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2661765) q[1];
sx q[1];
rz(-2.3485942) q[1];
sx q[1];
rz(-2.7686053) q[1];
rz(-pi) q[2];
rz(-0.21680253) q[3];
sx q[3];
rz(-0.6908145) q[3];
sx q[3];
rz(-1.3135214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5375157) q[2];
sx q[2];
rz(-2.5274369) q[2];
sx q[2];
rz(-2.3094731) q[2];
rz(2.4845947) q[3];
sx q[3];
rz(-1.4703581) q[3];
sx q[3];
rz(2.7895797) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.326139) q[0];
sx q[0];
rz(-0.6093381) q[0];
sx q[0];
rz(-2.4916008) q[0];
rz(0.12282148) q[1];
sx q[1];
rz(-0.42418066) q[1];
sx q[1];
rz(2.4688683) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8433326) q[0];
sx q[0];
rz(-1.6125868) q[0];
sx q[0];
rz(0.01173794) q[0];
rz(-pi) q[1];
rz(0.71750516) q[2];
sx q[2];
rz(-2.7360672) q[2];
sx q[2];
rz(0.44446352) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6482268) q[1];
sx q[1];
rz(-1.0208482) q[1];
sx q[1];
rz(2.9112766) q[1];
rz(-1.0510212) q[3];
sx q[3];
rz(-1.5477763) q[3];
sx q[3];
rz(-0.94205233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1136721) q[2];
sx q[2];
rz(-1.547926) q[2];
sx q[2];
rz(0.38635722) q[2];
rz(-1.1682642) q[3];
sx q[3];
rz(-1.9100274) q[3];
sx q[3];
rz(2.0333576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89762178) q[0];
sx q[0];
rz(-1.6672927) q[0];
sx q[0];
rz(3.0834055) q[0];
rz(2.8367786) q[1];
sx q[1];
rz(-2.0084281) q[1];
sx q[1];
rz(-2.286639) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8182504) q[0];
sx q[0];
rz(-1.5691681) q[0];
sx q[0];
rz(0.068308612) q[0];
rz(-2.0570175) q[2];
sx q[2];
rz(-2.2006319) q[2];
sx q[2];
rz(1.4818918) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82423174) q[1];
sx q[1];
rz(-2.3488725) q[1];
sx q[1];
rz(0.2893682) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0157973) q[3];
sx q[3];
rz(-1.7966174) q[3];
sx q[3];
rz(2.8005536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.052224) q[2];
sx q[2];
rz(-1.8171909) q[2];
sx q[2];
rz(-1.268187) q[2];
rz(0.44102272) q[3];
sx q[3];
rz(-2.5206168) q[3];
sx q[3];
rz(-2.2997901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(0.12255254) q[0];
sx q[0];
rz(-0.95128107) q[0];
sx q[0];
rz(2.4265491) q[0];
rz(2.7252281) q[1];
sx q[1];
rz(-2.6519897) q[1];
sx q[1];
rz(-2.8474836) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9884777) q[0];
sx q[0];
rz(-1.5128969) q[0];
sx q[0];
rz(0.87429509) q[0];
x q[1];
rz(-1.7753698) q[2];
sx q[2];
rz(-1.6343717) q[2];
sx q[2];
rz(2.3109189) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8079559) q[1];
sx q[1];
rz(-0.2273493) q[1];
sx q[1];
rz(-0.61509404) q[1];
rz(-pi) q[2];
rz(0.99331345) q[3];
sx q[3];
rz(-1.0619831) q[3];
sx q[3];
rz(-1.6779382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2405582) q[2];
sx q[2];
rz(-1.3996539) q[2];
sx q[2];
rz(-0.72875363) q[2];
rz(0.48480222) q[3];
sx q[3];
rz(-2.1383643) q[3];
sx q[3];
rz(0.17759594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4575551) q[0];
sx q[0];
rz(-1.8192679) q[0];
sx q[0];
rz(0.27873248) q[0];
rz(0.067525603) q[1];
sx q[1];
rz(-2.4261256) q[1];
sx q[1];
rz(2.6356437) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26209575) q[0];
sx q[0];
rz(-2.580632) q[0];
sx q[0];
rz(0.13763388) q[0];
rz(-2.9437482) q[2];
sx q[2];
rz(-1.4711507) q[2];
sx q[2];
rz(1.0464335) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5927148) q[1];
sx q[1];
rz(-2.1134106) q[1];
sx q[1];
rz(-0.9701258) q[1];
rz(-pi) q[2];
rz(1.9002731) q[3];
sx q[3];
rz(-2.4412324) q[3];
sx q[3];
rz(2.2735689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.94023306) q[2];
sx q[2];
rz(-1.6910005) q[2];
sx q[2];
rz(0.10227164) q[2];
rz(-0.30397948) q[3];
sx q[3];
rz(-0.18054466) q[3];
sx q[3];
rz(-0.48039082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2893696) q[0];
sx q[0];
rz(-2.3134573) q[0];
sx q[0];
rz(2.4612259) q[0];
rz(-0.96087372) q[1];
sx q[1];
rz(-1.5545132) q[1];
sx q[1];
rz(2.5804677) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2645123) q[0];
sx q[0];
rz(-1.3662369) q[0];
sx q[0];
rz(1.2007942) q[0];
rz(-0.84662171) q[2];
sx q[2];
rz(-1.4145803) q[2];
sx q[2];
rz(0.59000096) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.42237709) q[1];
sx q[1];
rz(-1.869162) q[1];
sx q[1];
rz(1.6279312) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64294358) q[3];
sx q[3];
rz(-0.74952945) q[3];
sx q[3];
rz(-1.7484322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23201021) q[2];
sx q[2];
rz(-2.1640615) q[2];
sx q[2];
rz(1.5634465) q[2];
rz(-2.1252508) q[3];
sx q[3];
rz(-1.7137824) q[3];
sx q[3];
rz(-2.0411172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72961724) q[0];
sx q[0];
rz(-0.53891861) q[0];
sx q[0];
rz(-2.6475661) q[0];
rz(1.5771075) q[1];
sx q[1];
rz(-0.99948519) q[1];
sx q[1];
rz(3.0513501) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099797332) q[0];
sx q[0];
rz(-1.8043515) q[0];
sx q[0];
rz(-2.6982624) q[0];
x q[1];
rz(0.87553067) q[2];
sx q[2];
rz(-1.1644496) q[2];
sx q[2];
rz(1.4515431) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4524432) q[1];
sx q[1];
rz(-1.1691368) q[1];
sx q[1];
rz(-2.8287877) q[1];
rz(-1.8537515) q[3];
sx q[3];
rz(-1.3850901) q[3];
sx q[3];
rz(2.9880092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.51284853) q[2];
sx q[2];
rz(-2.5675842) q[2];
sx q[2];
rz(2.4981456) q[2];
rz(-2.511034) q[3];
sx q[3];
rz(-0.75391155) q[3];
sx q[3];
rz(-3.0923617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028320463) q[0];
sx q[0];
rz(-1.4406818) q[0];
sx q[0];
rz(1.2233618) q[0];
rz(-2.2471097) q[1];
sx q[1];
rz(-0.62791413) q[1];
sx q[1];
rz(1.1462513) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.664331) q[0];
sx q[0];
rz(-1.0205752) q[0];
sx q[0];
rz(-2.6653637) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9983534) q[2];
sx q[2];
rz(-1.9263066) q[2];
sx q[2];
rz(1.8380788) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8889956) q[1];
sx q[1];
rz(-1.3301783) q[1];
sx q[1];
rz(0.29901286) q[1];
x q[2];
rz(0.98553879) q[3];
sx q[3];
rz(-1.7904591) q[3];
sx q[3];
rz(-2.1287763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2724096) q[2];
sx q[2];
rz(-0.86239186) q[2];
sx q[2];
rz(1.4746846) q[2];
rz(0.21371755) q[3];
sx q[3];
rz(-1.7361879) q[3];
sx q[3];
rz(0.86764446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2992582) q[0];
sx q[0];
rz(-2.5304351) q[0];
sx q[0];
rz(-2.4364731) q[0];
rz(-0.57535386) q[1];
sx q[1];
rz(-1.7333142) q[1];
sx q[1];
rz(-2.5720678) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0234719) q[0];
sx q[0];
rz(-3.0295353) q[0];
sx q[0];
rz(2.3311291) q[0];
x q[1];
rz(0.51409431) q[2];
sx q[2];
rz(-0.55214685) q[2];
sx q[2];
rz(-0.29468003) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8882506) q[1];
sx q[1];
rz(-0.98501316) q[1];
sx q[1];
rz(-0.82510494) q[1];
rz(-pi) q[2];
x q[2];
rz(0.488987) q[3];
sx q[3];
rz(-1.8400116) q[3];
sx q[3];
rz(-1.0753635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.19499714) q[2];
sx q[2];
rz(-0.94299287) q[2];
sx q[2];
rz(0.97879624) q[2];
rz(0.43092522) q[3];
sx q[3];
rz(-0.48723358) q[3];
sx q[3];
rz(1.1270969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2836714) q[0];
sx q[0];
rz(-3.1223174) q[0];
sx q[0];
rz(1.4625782) q[0];
rz(-1.0390394) q[1];
sx q[1];
rz(-1.3835399) q[1];
sx q[1];
rz(-1.9336112) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9383065) q[0];
sx q[0];
rz(-1.1933367) q[0];
sx q[0];
rz(-2.8164144) q[0];
rz(-pi) q[1];
rz(-2.0956089) q[2];
sx q[2];
rz(-1.599858) q[2];
sx q[2];
rz(-0.20653221) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4962311) q[1];
sx q[1];
rz(-1.4252311) q[1];
sx q[1];
rz(0.95349378) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7321587) q[3];
sx q[3];
rz(-0.40009016) q[3];
sx q[3];
rz(2.0619947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.98564467) q[2];
sx q[2];
rz(-1.0040823) q[2];
sx q[2];
rz(-1.3204302) q[2];
rz(-2.7909347) q[3];
sx q[3];
rz(-0.81736332) q[3];
sx q[3];
rz(-0.58026522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88631267) q[0];
sx q[0];
rz(-1.9552312) q[0];
sx q[0];
rz(2.2807688) q[0];
rz(1.6571922) q[1];
sx q[1];
rz(-1.7488372) q[1];
sx q[1];
rz(-1.7120672) q[1];
rz(-1.4996573) q[2];
sx q[2];
rz(-1.9491458) q[2];
sx q[2];
rz(-0.083935621) q[2];
rz(-2.8410925) q[3];
sx q[3];
rz(-0.90682744) q[3];
sx q[3];
rz(-1.711267) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
