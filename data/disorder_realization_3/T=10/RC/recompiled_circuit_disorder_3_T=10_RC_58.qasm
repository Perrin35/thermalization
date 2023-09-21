OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.91575032) q[0];
sx q[0];
rz(3.1728035) q[0];
sx q[0];
rz(6.7682545) q[0];
rz(0.78753161) q[1];
sx q[1];
rz(-1.0163611) q[1];
sx q[1];
rz(2.7273942) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6611377) q[0];
sx q[0];
rz(-1.2167131) q[0];
sx q[0];
rz(0.35868355) q[0];
rz(-3.0468416) q[2];
sx q[2];
rz(-2.13846) q[2];
sx q[2];
rz(1.3324347) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1297076) q[1];
sx q[1];
rz(-1.9958479) q[1];
sx q[1];
rz(-3.0711864) q[1];
x q[2];
rz(1.6730509) q[3];
sx q[3];
rz(-1.8342606) q[3];
sx q[3];
rz(1.9139569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2177314) q[2];
sx q[2];
rz(-1.8863181) q[2];
sx q[2];
rz(-3.1100173) q[2];
rz(-1.2565553) q[3];
sx q[3];
rz(-2.6387408) q[3];
sx q[3];
rz(2.6149635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.8388222) q[0];
sx q[0];
rz(-1.6844203) q[0];
sx q[0];
rz(-2.9717428) q[0];
rz(-0.70392144) q[1];
sx q[1];
rz(-2.070065) q[1];
sx q[1];
rz(-2.6020715) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3300433) q[0];
sx q[0];
rz(-1.1216315) q[0];
sx q[0];
rz(3.0898068) q[0];
x q[1];
rz(1.1621446) q[2];
sx q[2];
rz(-0.89102972) q[2];
sx q[2];
rz(1.0811999) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.62944618) q[1];
sx q[1];
rz(-2.0680032) q[1];
sx q[1];
rz(-2.9260103) q[1];
rz(2.5856421) q[3];
sx q[3];
rz(-2.2009938) q[3];
sx q[3];
rz(-2.6549908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4743621) q[2];
sx q[2];
rz(-1.2381866) q[2];
sx q[2];
rz(-2.898522) q[2];
rz(0.66611755) q[3];
sx q[3];
rz(-0.56454286) q[3];
sx q[3];
rz(1.8977785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3617525) q[0];
sx q[0];
rz(-0.11479522) q[0];
sx q[0];
rz(-2.6932122) q[0];
rz(-1.7547296) q[1];
sx q[1];
rz(-1.9877537) q[1];
sx q[1];
rz(0.2562491) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019779531) q[0];
sx q[0];
rz(-2.461713) q[0];
sx q[0];
rz(-2.7133184) q[0];
x q[1];
rz(2.0492378) q[2];
sx q[2];
rz(-1.97176) q[2];
sx q[2];
rz(-1.5329597) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5970522) q[1];
sx q[1];
rz(-1.2433194) q[1];
sx q[1];
rz(-2.271133) q[1];
rz(2.933421) q[3];
sx q[3];
rz(-2.9609207) q[3];
sx q[3];
rz(0.65073035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3391352) q[2];
sx q[2];
rz(-0.70636237) q[2];
sx q[2];
rz(0.57470542) q[2];
rz(-1.7859219) q[3];
sx q[3];
rz(-1.1698497) q[3];
sx q[3];
rz(-1.908196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9280076) q[0];
sx q[0];
rz(-1.428823) q[0];
sx q[0];
rz(2.8821049) q[0];
rz(-1.9909987) q[1];
sx q[1];
rz(-1.3535627) q[1];
sx q[1];
rz(-2.4096699) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5314732) q[0];
sx q[0];
rz(-1.1024794) q[0];
sx q[0];
rz(-2.5584695) q[0];
x q[1];
rz(-0.133693) q[2];
sx q[2];
rz(-0.72172726) q[2];
sx q[2];
rz(-2.0267817) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.54938984) q[1];
sx q[1];
rz(-1.7648186) q[1];
sx q[1];
rz(0.26280304) q[1];
rz(-pi) q[2];
rz(-0.22051208) q[3];
sx q[3];
rz(-1.4482499) q[3];
sx q[3];
rz(-1.8054655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6966454) q[2];
sx q[2];
rz(-1.8680633) q[2];
sx q[2];
rz(-0.15110061) q[2];
rz(0.54667306) q[3];
sx q[3];
rz(-1.0497382) q[3];
sx q[3];
rz(0.37724075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-0.54995173) q[0];
sx q[0];
rz(-1.0629835) q[0];
sx q[0];
rz(-1.8792101) q[0];
rz(1.6732015) q[1];
sx q[1];
rz(-2.5322798) q[1];
sx q[1];
rz(-2.343822) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20194963) q[0];
sx q[0];
rz(-1.573283) q[0];
sx q[0];
rz(-1.6909088) q[0];
rz(-pi) q[1];
rz(-0.61879976) q[2];
sx q[2];
rz(-0.36703645) q[2];
sx q[2];
rz(2.7235081) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5835727) q[1];
sx q[1];
rz(-2.1961374) q[1];
sx q[1];
rz(0.11548345) q[1];
rz(2.5161414) q[3];
sx q[3];
rz(-2.6126385) q[3];
sx q[3];
rz(-1.7367425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.795934) q[2];
sx q[2];
rz(-2.5107333) q[2];
sx q[2];
rz(-2.8395555) q[2];
rz(1.1473514) q[3];
sx q[3];
rz(-1.6796422) q[3];
sx q[3];
rz(-0.54774493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40925947) q[0];
sx q[0];
rz(-3.0497666) q[0];
sx q[0];
rz(-1.1557895) q[0];
rz(-1.0844768) q[1];
sx q[1];
rz(-2.1613354) q[1];
sx q[1];
rz(-0.070080431) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87175831) q[0];
sx q[0];
rz(-1.4176798) q[0];
sx q[0];
rz(-1.7437177) q[0];
rz(-2.5885133) q[2];
sx q[2];
rz(-0.86420176) q[2];
sx q[2];
rz(1.4097708) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3190223) q[1];
sx q[1];
rz(-0.59148568) q[1];
sx q[1];
rz(-0.2925847) q[1];
rz(1.6453679) q[3];
sx q[3];
rz(-1.5354904) q[3];
sx q[3];
rz(1.1535742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8391116) q[2];
sx q[2];
rz(-1.9589067) q[2];
sx q[2];
rz(2.5202259) q[2];
rz(1.4403884) q[3];
sx q[3];
rz(-2.6337603) q[3];
sx q[3];
rz(0.5733718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0550585) q[0];
sx q[0];
rz(-1.4165514) q[0];
sx q[0];
rz(2.281718) q[0];
rz(1.9372008) q[1];
sx q[1];
rz(-2.269373) q[1];
sx q[1];
rz(0.0079356114) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.313109) q[0];
sx q[0];
rz(-2.5223753) q[0];
sx q[0];
rz(-0.45066582) q[0];
rz(-pi) q[1];
rz(-0.87848778) q[2];
sx q[2];
rz(-1.9108859) q[2];
sx q[2];
rz(-2.0932587) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7759526) q[1];
sx q[1];
rz(-2.0740293) q[1];
sx q[1];
rz(-2.3535437) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2357335) q[3];
sx q[3];
rz(-2.6568036) q[3];
sx q[3];
rz(-1.0672027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.69592151) q[2];
sx q[2];
rz(-1.9182703) q[2];
sx q[2];
rz(-2.725214) q[2];
rz(-1.773206) q[3];
sx q[3];
rz(-1.2973283) q[3];
sx q[3];
rz(-0.90464512) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96173441) q[0];
sx q[0];
rz(-0.27357736) q[0];
sx q[0];
rz(0.36488786) q[0];
rz(2.2015613) q[1];
sx q[1];
rz(-0.53661984) q[1];
sx q[1];
rz(1.6392802) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84425981) q[0];
sx q[0];
rz(-1.5607921) q[0];
sx q[0];
rz(2.8914333) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47183581) q[2];
sx q[2];
rz(-0.89315692) q[2];
sx q[2];
rz(0.0080646947) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.18174905) q[1];
sx q[1];
rz(-1.0105003) q[1];
sx q[1];
rz(-0.77397857) q[1];
rz(-pi) q[2];
rz(1.3650465) q[3];
sx q[3];
rz(-1.4806517) q[3];
sx q[3];
rz(2.3004325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.49729785) q[2];
sx q[2];
rz(-2.6361894) q[2];
sx q[2];
rz(1.0650744) q[2];
rz(-0.30125695) q[3];
sx q[3];
rz(-1.4091636) q[3];
sx q[3];
rz(-1.6206954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.168468) q[0];
sx q[0];
rz(-3.0597866) q[0];
sx q[0];
rz(-2.705943) q[0];
rz(-1.7565953) q[1];
sx q[1];
rz(-0.47043097) q[1];
sx q[1];
rz(0.41697821) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5676253) q[0];
sx q[0];
rz(-1.6084533) q[0];
sx q[0];
rz(2.2249939) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66395335) q[2];
sx q[2];
rz(-1.10154) q[2];
sx q[2];
rz(0.19367733) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7855362) q[1];
sx q[1];
rz(-1.0724663) q[1];
sx q[1];
rz(-0.15798012) q[1];
rz(-pi) q[2];
rz(0.74420332) q[3];
sx q[3];
rz(-2.9508698) q[3];
sx q[3];
rz(-0.71803367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.10432648) q[2];
sx q[2];
rz(-1.4929079) q[2];
sx q[2];
rz(-2.9157675) q[2];
rz(-2.9337692) q[3];
sx q[3];
rz(-0.72312975) q[3];
sx q[3];
rz(-2.5549755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.726783) q[0];
sx q[0];
rz(-0.2668969) q[0];
sx q[0];
rz(1.6171932) q[0];
rz(-2.1879451) q[1];
sx q[1];
rz(-1.2748268) q[1];
sx q[1];
rz(1.8189925) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7693217) q[0];
sx q[0];
rz(-1.9479381) q[0];
sx q[0];
rz(-3.0110714) q[0];
x q[1];
rz(-2.0595423) q[2];
sx q[2];
rz(-0.49968038) q[2];
sx q[2];
rz(0.79364712) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0304071) q[1];
sx q[1];
rz(-0.92799312) q[1];
sx q[1];
rz(-2.2305957) q[1];
rz(-pi) q[2];
rz(0.91536509) q[3];
sx q[3];
rz(-1.4487106) q[3];
sx q[3];
rz(-1.1673499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7913197) q[2];
sx q[2];
rz(-2.0952756) q[2];
sx q[2];
rz(-1.0160758) q[2];
rz(1.919205) q[3];
sx q[3];
rz(-0.17761579) q[3];
sx q[3];
rz(-2.5861752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9983457) q[0];
sx q[0];
rz(-0.9220985) q[0];
sx q[0];
rz(-2.0621598) q[0];
rz(1.3636419) q[1];
sx q[1];
rz(-1.9121871) q[1];
sx q[1];
rz(-1.8180064) q[1];
rz(-0.50921847) q[2];
sx q[2];
rz(-1.5621395) q[2];
sx q[2];
rz(-0.12315673) q[2];
rz(-1.3397459) q[3];
sx q[3];
rz(-1.6587202) q[3];
sx q[3];
rz(-0.29826577) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];