OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0619573) q[0];
sx q[0];
rz(-0.24793967) q[0];
sx q[0];
rz(2.3214582) q[0];
rz(-3.1007383) q[1];
sx q[1];
rz(-0.78894579) q[1];
sx q[1];
rz(-0.087892428) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0996617) q[0];
sx q[0];
rz(-2.4798205) q[0];
sx q[0];
rz(-2.3753138) q[0];
rz(0.90859969) q[2];
sx q[2];
rz(-0.63630051) q[2];
sx q[2];
rz(1.5719906) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3327515) q[1];
sx q[1];
rz(-1.0665227) q[1];
sx q[1];
rz(1.5682975) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8294931) q[3];
sx q[3];
rz(-1.5651349) q[3];
sx q[3];
rz(2.2744846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1408046) q[2];
sx q[2];
rz(-1.452383) q[2];
sx q[2];
rz(-0.56935707) q[2];
rz(1.6128444) q[3];
sx q[3];
rz(-0.6362392) q[3];
sx q[3];
rz(1.7830085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59812087) q[0];
sx q[0];
rz(-0.16508979) q[0];
sx q[0];
rz(-0.55364451) q[0];
rz(1.9042632) q[1];
sx q[1];
rz(-1.7747223) q[1];
sx q[1];
rz(-1.233261) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.935826) q[0];
sx q[0];
rz(-2.2201949) q[0];
sx q[0];
rz(-0.68769023) q[0];
rz(2.4670062) q[2];
sx q[2];
rz(-1.6709176) q[2];
sx q[2];
rz(1.1652201) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1439647) q[1];
sx q[1];
rz(-1.7768754) q[1];
sx q[1];
rz(-3.1254083) q[1];
x q[2];
rz(-2.8418856) q[3];
sx q[3];
rz(-2.4518659) q[3];
sx q[3];
rz(1.0563861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1374986) q[2];
sx q[2];
rz(-1.4592905) q[2];
sx q[2];
rz(-0.0022350524) q[2];
rz(0.8301174) q[3];
sx q[3];
rz(-2.3486962) q[3];
sx q[3];
rz(1.2494276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24901351) q[0];
sx q[0];
rz(-0.61616388) q[0];
sx q[0];
rz(2.2900443) q[0];
rz(2.3705204) q[1];
sx q[1];
rz(-1.4665736) q[1];
sx q[1];
rz(3.070389) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035633798) q[0];
sx q[0];
rz(-1.7674315) q[0];
sx q[0];
rz(1.8133624) q[0];
x q[1];
rz(-2.0257646) q[2];
sx q[2];
rz(-1.8075382) q[2];
sx q[2];
rz(0.67982212) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5350266) q[1];
sx q[1];
rz(-2.2177794) q[1];
sx q[1];
rz(1.4024629) q[1];
rz(-0.72910587) q[3];
sx q[3];
rz(-0.90568542) q[3];
sx q[3];
rz(0.8641181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3123793) q[2];
sx q[2];
rz(-0.92386121) q[2];
sx q[2];
rz(-0.52345792) q[2];
rz(-0.33189014) q[3];
sx q[3];
rz(-0.060398014) q[3];
sx q[3];
rz(0.50326842) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7317384) q[0];
sx q[0];
rz(-2.8926352) q[0];
sx q[0];
rz(-0.82558924) q[0];
rz(1.9748953) q[1];
sx q[1];
rz(-0.86667934) q[1];
sx q[1];
rz(1.694214) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4226845) q[0];
sx q[0];
rz(-1.3410765) q[0];
sx q[0];
rz(-2.9740888) q[0];
rz(-pi) q[1];
rz(-1.3924613) q[2];
sx q[2];
rz(-1.5797838) q[2];
sx q[2];
rz(-2.044878) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4923258) q[1];
sx q[1];
rz(-2.4267303) q[1];
sx q[1];
rz(1.3267172) q[1];
rz(-2.7778266) q[3];
sx q[3];
rz(-1.3437265) q[3];
sx q[3];
rz(2.098339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6772785) q[2];
sx q[2];
rz(-1.7769287) q[2];
sx q[2];
rz(-2.9966667) q[2];
rz(-2.1302917) q[3];
sx q[3];
rz(-0.62711182) q[3];
sx q[3];
rz(-3.1269126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1458364) q[0];
sx q[0];
rz(-0.053481426) q[0];
sx q[0];
rz(2.4575535) q[0];
rz(1.1902635) q[1];
sx q[1];
rz(-0.66527706) q[1];
sx q[1];
rz(-1.9285944) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1367462) q[0];
sx q[0];
rz(-0.46127013) q[0];
sx q[0];
rz(-0.61704163) q[0];
rz(2.331779) q[2];
sx q[2];
rz(-2.3286616) q[2];
sx q[2];
rz(-0.44250689) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.1274384) q[1];
sx q[1];
rz(-1.6886097) q[1];
sx q[1];
rz(0.25728667) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1608859) q[3];
sx q[3];
rz(-0.31168918) q[3];
sx q[3];
rz(-2.2463472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6325536) q[2];
sx q[2];
rz(-1.2712487) q[2];
sx q[2];
rz(-0.63009134) q[2];
rz(0.57224327) q[3];
sx q[3];
rz(-2.4961491) q[3];
sx q[3];
rz(0.8852638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0260139) q[0];
sx q[0];
rz(-2.2981839) q[0];
sx q[0];
rz(-0.28636006) q[0];
rz(-0.59965602) q[1];
sx q[1];
rz(-1.1279794) q[1];
sx q[1];
rz(-2.5568331) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.023914) q[0];
sx q[0];
rz(-0.6891793) q[0];
sx q[0];
rz(0.50424772) q[0];
x q[1];
rz(1.5445968) q[2];
sx q[2];
rz(-0.95653406) q[2];
sx q[2];
rz(-2.7363077) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.11440052) q[1];
sx q[1];
rz(-0.90585867) q[1];
sx q[1];
rz(0.36892885) q[1];
rz(2.9752475) q[3];
sx q[3];
rz(-1.5208941) q[3];
sx q[3];
rz(2.9970616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9770603) q[2];
sx q[2];
rz(-2.2571199) q[2];
sx q[2];
rz(-1.6284846) q[2];
rz(0.96308723) q[3];
sx q[3];
rz(-1.3685127) q[3];
sx q[3];
rz(2.517038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088783711) q[0];
sx q[0];
rz(-1.7099021) q[0];
sx q[0];
rz(-0.54876304) q[0];
rz(-0.051963003) q[1];
sx q[1];
rz(-0.49182645) q[1];
sx q[1];
rz(2.5351977) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43019766) q[0];
sx q[0];
rz(-2.093962) q[0];
sx q[0];
rz(2.8682312) q[0];
x q[1];
rz(2.0386001) q[2];
sx q[2];
rz(-0.71873795) q[2];
sx q[2];
rz(1.7407835) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.079271) q[1];
sx q[1];
rz(-0.81264001) q[1];
sx q[1];
rz(2.9864242) q[1];
x q[2];
rz(2.7860836) q[3];
sx q[3];
rz(-1.9572557) q[3];
sx q[3];
rz(-0.10458065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9397883) q[2];
sx q[2];
rz(-0.99210343) q[2];
sx q[2];
rz(-2.6331804) q[2];
rz(-1.1172179) q[3];
sx q[3];
rz(-2.9682187) q[3];
sx q[3];
rz(1.6569051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8023119) q[0];
sx q[0];
rz(-2.7957714) q[0];
sx q[0];
rz(-0.75396496) q[0];
rz(2.1915961) q[1];
sx q[1];
rz(-1.6597304) q[1];
sx q[1];
rz(-2.9139013) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75922155) q[0];
sx q[0];
rz(-1.641944) q[0];
sx q[0];
rz(1.9842149) q[0];
rz(-pi) q[1];
rz(1.9367427) q[2];
sx q[2];
rz(-1.9908675) q[2];
sx q[2];
rz(2.7811188) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5880809) q[1];
sx q[1];
rz(-1.7953824) q[1];
sx q[1];
rz(-2.5976546) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5826107) q[3];
sx q[3];
rz(-1.535927) q[3];
sx q[3];
rz(0.99136409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0630539) q[2];
sx q[2];
rz(-0.20919122) q[2];
sx q[2];
rz(0.83855808) q[2];
rz(2.7770384) q[3];
sx q[3];
rz(-1.3804599) q[3];
sx q[3];
rz(2.3013039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7446328) q[0];
sx q[0];
rz(-1.3752221) q[0];
sx q[0];
rz(2.3378085) q[0];
rz(1.0546168) q[1];
sx q[1];
rz(-2.5024253) q[1];
sx q[1];
rz(-1.1484336) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6156857) q[0];
sx q[0];
rz(-1.4988168) q[0];
sx q[0];
rz(3.1176438) q[0];
rz(2.7788413) q[2];
sx q[2];
rz(-2.0667549) q[2];
sx q[2];
rz(-2.9438058) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2003186) q[1];
sx q[1];
rz(-2.4793844) q[1];
sx q[1];
rz(0.83076417) q[1];
rz(-2.7896342) q[3];
sx q[3];
rz(-1.6037914) q[3];
sx q[3];
rz(-1.2034504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.049872963) q[2];
sx q[2];
rz(-1.0174948) q[2];
sx q[2];
rz(2.3633374) q[2];
rz(-0.19566472) q[3];
sx q[3];
rz(-1.0396495) q[3];
sx q[3];
rz(-1.0218609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2739928) q[0];
sx q[0];
rz(-0.99853981) q[0];
sx q[0];
rz(-2.8883949) q[0];
rz(2.4018438) q[1];
sx q[1];
rz(-2.4052129) q[1];
sx q[1];
rz(0.69828066) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2525576) q[0];
sx q[0];
rz(-1.4989984) q[0];
sx q[0];
rz(-0.30307146) q[0];
rz(1.8928705) q[2];
sx q[2];
rz(-0.078157166) q[2];
sx q[2];
rz(-0.31101481) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.89430289) q[1];
sx q[1];
rz(-1.3548046) q[1];
sx q[1];
rz(1.1029907) q[1];
rz(1.2020338) q[3];
sx q[3];
rz(-2.5420815) q[3];
sx q[3];
rz(-2.5332019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.11761052) q[2];
sx q[2];
rz(-0.88279804) q[2];
sx q[2];
rz(-0.84021604) q[2];
rz(-0.64030567) q[3];
sx q[3];
rz(-0.87602031) q[3];
sx q[3];
rz(1.1673814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2587851) q[0];
sx q[0];
rz(-2.2866645) q[0];
sx q[0];
rz(-1.4186161) q[0];
rz(-3.1124658) q[1];
sx q[1];
rz(-0.11321414) q[1];
sx q[1];
rz(-1.3197457) q[1];
rz(1.9188948) q[2];
sx q[2];
rz(-2.1445027) q[2];
sx q[2];
rz(2.4377844) q[2];
rz(-0.15016951) q[3];
sx q[3];
rz(-0.88610813) q[3];
sx q[3];
rz(-0.025972493) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];