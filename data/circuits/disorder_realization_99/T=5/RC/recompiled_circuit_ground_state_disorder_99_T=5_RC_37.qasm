OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.91166624) q[0];
sx q[0];
rz(-0.32719964) q[0];
sx q[0];
rz(1.3258452) q[0];
rz(-1.7757379) q[1];
sx q[1];
rz(-2.7471625) q[1];
sx q[1];
rz(2.3529513) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45730725) q[0];
sx q[0];
rz(-2.4198902) q[0];
sx q[0];
rz(-1.0941605) q[0];
rz(-pi) q[1];
rz(2.2872988) q[2];
sx q[2];
rz(-2.2546356) q[2];
sx q[2];
rz(1.270806) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0683191) q[1];
sx q[1];
rz(-1.8536659) q[1];
sx q[1];
rz(-1.7712996) q[1];
x q[2];
rz(1.7472505) q[3];
sx q[3];
rz(-0.6354699) q[3];
sx q[3];
rz(0.42144767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6476562) q[2];
sx q[2];
rz(-1.4467354) q[2];
sx q[2];
rz(2.9392865) q[2];
rz(3.0322187) q[3];
sx q[3];
rz(-2.5325363) q[3];
sx q[3];
rz(3.0561395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.0789455) q[0];
sx q[0];
rz(-2.2181692) q[0];
sx q[0];
rz(-1.1751291) q[0];
rz(-1.6773978) q[1];
sx q[1];
rz(-1.0025832) q[1];
sx q[1];
rz(-0.35951231) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8749411) q[0];
sx q[0];
rz(-1.4066937) q[0];
sx q[0];
rz(0.070673857) q[0];
rz(0.56698842) q[2];
sx q[2];
rz(-2.5051281) q[2];
sx q[2];
rz(1.3186962) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2613758) q[1];
sx q[1];
rz(-1.2296074) q[1];
sx q[1];
rz(0.64067322) q[1];
x q[2];
rz(-0.58200804) q[3];
sx q[3];
rz(-1.0543543) q[3];
sx q[3];
rz(-0.88508115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.29490617) q[2];
sx q[2];
rz(-1.0470231) q[2];
sx q[2];
rz(-2.8666551) q[2];
rz(-1.8168195) q[3];
sx q[3];
rz(-1.5271657) q[3];
sx q[3];
rz(-1.2980609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0637958) q[0];
sx q[0];
rz(-0.15811385) q[0];
sx q[0];
rz(-2.7445444) q[0];
rz(-0.60931698) q[1];
sx q[1];
rz(-1.3343697) q[1];
sx q[1];
rz(-1.5416001) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5865094) q[0];
sx q[0];
rz(-1.6753046) q[0];
sx q[0];
rz(-1.2782607) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3808525) q[2];
sx q[2];
rz(-0.06643387) q[2];
sx q[2];
rz(-2.9125467) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0761281) q[1];
sx q[1];
rz(-2.4148126) q[1];
sx q[1];
rz(-1.5453592) q[1];
rz(-2.479781) q[3];
sx q[3];
rz(-0.99259171) q[3];
sx q[3];
rz(-0.44826642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2043173) q[2];
sx q[2];
rz(-1.5466362) q[2];
sx q[2];
rz(-0.23359648) q[2];
rz(-1.8448081) q[3];
sx q[3];
rz(-1.9498884) q[3];
sx q[3];
rz(-0.72062033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63491708) q[0];
sx q[0];
rz(-1.9838061) q[0];
sx q[0];
rz(-1.3651715) q[0];
rz(-1.8695976) q[1];
sx q[1];
rz(-1.9422266) q[1];
sx q[1];
rz(-1.2619527) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8376217) q[0];
sx q[0];
rz(-0.6084992) q[0];
sx q[0];
rz(-2.801218) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6532516) q[2];
sx q[2];
rz(-2.8318498) q[2];
sx q[2];
rz(2.7249634) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8820928) q[1];
sx q[1];
rz(-1.6354523) q[1];
sx q[1];
rz(0.18292173) q[1];
rz(0.19660321) q[3];
sx q[3];
rz(-2.3075309) q[3];
sx q[3];
rz(-1.3722668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9103553) q[2];
sx q[2];
rz(-2.1114025) q[2];
sx q[2];
rz(-2.3528698) q[2];
rz(-1.8606868) q[3];
sx q[3];
rz(-1.3642045) q[3];
sx q[3];
rz(-2.1196712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(1.7541517) q[0];
sx q[0];
rz(-1.0197637) q[0];
sx q[0];
rz(-2.4805241) q[0];
rz(-2.0263653) q[1];
sx q[1];
rz(-1.8293569) q[1];
sx q[1];
rz(0.95975319) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57908995) q[0];
sx q[0];
rz(-1.8402365) q[0];
sx q[0];
rz(0.53621063) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16483817) q[2];
sx q[2];
rz(-1.0292813) q[2];
sx q[2];
rz(-2.7680226) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.58127357) q[1];
sx q[1];
rz(-1.7035824) q[1];
sx q[1];
rz(1.2788692) q[1];
rz(0.19790217) q[3];
sx q[3];
rz(-0.93557916) q[3];
sx q[3];
rz(2.9445348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.58012086) q[2];
sx q[2];
rz(-2.6802907) q[2];
sx q[2];
rz(-2.0816154) q[2];
rz(-2.3914242) q[3];
sx q[3];
rz(-1.7629905) q[3];
sx q[3];
rz(-1.9157971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8691413) q[0];
sx q[0];
rz(-1.5831818) q[0];
sx q[0];
rz(0.22431746) q[0];
rz(1.6250826) q[1];
sx q[1];
rz(-2.5254011) q[1];
sx q[1];
rz(-3.065965) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.513994) q[0];
sx q[0];
rz(-0.69116163) q[0];
sx q[0];
rz(-0.94011728) q[0];
rz(-pi) q[1];
rz(0.096616726) q[2];
sx q[2];
rz(-1.268309) q[2];
sx q[2];
rz(2.5563478) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5995248) q[1];
sx q[1];
rz(-1.1097399) q[1];
sx q[1];
rz(0.49054029) q[1];
rz(-pi) q[2];
rz(-0.43419713) q[3];
sx q[3];
rz(-0.89354529) q[3];
sx q[3];
rz(-1.7627258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.533796) q[2];
sx q[2];
rz(-0.7527315) q[2];
sx q[2];
rz(0.047018615) q[2];
rz(-2.4646344) q[3];
sx q[3];
rz(-1.2983863) q[3];
sx q[3];
rz(3.1162139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9852279) q[0];
sx q[0];
rz(-1.4893091) q[0];
sx q[0];
rz(1.4439247) q[0];
rz(-2.2185183) q[1];
sx q[1];
rz(-1.0052899) q[1];
sx q[1];
rz(-2.9583171) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5248337) q[0];
sx q[0];
rz(-1.3206894) q[0];
sx q[0];
rz(-0.5768187) q[0];
rz(-pi) q[1];
rz(-0.44281339) q[2];
sx q[2];
rz(-1.1842791) q[2];
sx q[2];
rz(0.3511951) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0456775) q[1];
sx q[1];
rz(-2.5325477) q[1];
sx q[1];
rz(0.97859971) q[1];
rz(-pi) q[2];
rz(0.28740643) q[3];
sx q[3];
rz(-2.1309521) q[3];
sx q[3];
rz(1.0792102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.034059374) q[2];
sx q[2];
rz(-1.6851765) q[2];
sx q[2];
rz(-3.122186) q[2];
rz(2.9803045) q[3];
sx q[3];
rz(-1.015377) q[3];
sx q[3];
rz(1.9136782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0374544) q[0];
sx q[0];
rz(-0.4087953) q[0];
sx q[0];
rz(3.0492875) q[0];
rz(-1.1760938) q[1];
sx q[1];
rz(-2.8832925) q[1];
sx q[1];
rz(1.4091122) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8831439) q[0];
sx q[0];
rz(-1.6250659) q[0];
sx q[0];
rz(1.6759713) q[0];
x q[1];
rz(-0.95887948) q[2];
sx q[2];
rz(-1.9643219) q[2];
sx q[2];
rz(0.78262353) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2829874) q[1];
sx q[1];
rz(-1.8388565) q[1];
sx q[1];
rz(0.40045935) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61586611) q[3];
sx q[3];
rz(-2.0558254) q[3];
sx q[3];
rz(-0.46014338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1000503) q[2];
sx q[2];
rz(-1.5771733) q[2];
sx q[2];
rz(1.9680295) q[2];
rz(-0.38343492) q[3];
sx q[3];
rz(-1.3078559) q[3];
sx q[3];
rz(-2.0508544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1293056) q[0];
sx q[0];
rz(-1.4947991) q[0];
sx q[0];
rz(-0.20558414) q[0];
rz(0.81941098) q[1];
sx q[1];
rz(-1.2148379) q[1];
sx q[1];
rz(-0.49088556) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6554479) q[0];
sx q[0];
rz(-1.1974338) q[0];
sx q[0];
rz(-2.4400737) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41342469) q[2];
sx q[2];
rz(-0.67081645) q[2];
sx q[2];
rz(2.6476423) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6632884) q[1];
sx q[1];
rz(-2.5712447) q[1];
sx q[1];
rz(1.4823518) q[1];
rz(-3.0428995) q[3];
sx q[3];
rz(-1.0617563) q[3];
sx q[3];
rz(-0.72580298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4129591) q[2];
sx q[2];
rz(-2.1042991) q[2];
sx q[2];
rz(1.884985) q[2];
rz(2.7058153) q[3];
sx q[3];
rz(-1.1083009) q[3];
sx q[3];
rz(-0.88424879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7213781) q[0];
sx q[0];
rz(-0.90737897) q[0];
sx q[0];
rz(-0.23289982) q[0];
rz(-0.13889343) q[1];
sx q[1];
rz(-1.1537617) q[1];
sx q[1];
rz(-0.5084261) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.229436) q[0];
sx q[0];
rz(-0.60750738) q[0];
sx q[0];
rz(-1.1077393) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0632238) q[2];
sx q[2];
rz(-1.3295577) q[2];
sx q[2];
rz(-2.9335748) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0838937) q[1];
sx q[1];
rz(-3.0073799) q[1];
sx q[1];
rz(-1.6300409) q[1];
rz(-pi) q[2];
rz(-2.2670161) q[3];
sx q[3];
rz(-1.1406058) q[3];
sx q[3];
rz(-0.065635292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7306708) q[2];
sx q[2];
rz(-1.4620917) q[2];
sx q[2];
rz(0.81653583) q[2];
rz(0.38980347) q[3];
sx q[3];
rz(-1.0023508) q[3];
sx q[3];
rz(1.3780814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.447406) q[0];
sx q[0];
rz(-0.92300713) q[0];
sx q[0];
rz(-0.22966455) q[0];
rz(1.6387088) q[1];
sx q[1];
rz(-1.3353744) q[1];
sx q[1];
rz(-2.2257805) q[1];
rz(-1.3982795) q[2];
sx q[2];
rz(-1.1926706) q[2];
sx q[2];
rz(-0.43422912) q[2];
rz(-2.2727454) q[3];
sx q[3];
rz(-1.6985536) q[3];
sx q[3];
rz(2.539195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
