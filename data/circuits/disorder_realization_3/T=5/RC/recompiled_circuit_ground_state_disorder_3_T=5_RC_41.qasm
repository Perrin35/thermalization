OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.42263) q[0];
sx q[0];
rz(-2.842272) q[0];
sx q[0];
rz(2.646995) q[0];
rz(-1.9994796) q[1];
sx q[1];
rz(-2.1358868) q[1];
sx q[1];
rz(-1.1297273) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.513213) q[0];
sx q[0];
rz(-2.2404379) q[0];
sx q[0];
rz(-1.7731401) q[0];
x q[1];
rz(2.3006265) q[2];
sx q[2];
rz(-0.6305002) q[2];
sx q[2];
rz(1.7918685) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2430302) q[1];
sx q[1];
rz(-2.3645325) q[1];
sx q[1];
rz(-1.9352566) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7143728) q[3];
sx q[3];
rz(-0.96517206) q[3];
sx q[3];
rz(2.1744414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8226681) q[2];
sx q[2];
rz(-1.319004) q[2];
sx q[2];
rz(2.2887716) q[2];
rz(-1.3301814) q[3];
sx q[3];
rz(-0.69555247) q[3];
sx q[3];
rz(0.046796355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80832076) q[0];
sx q[0];
rz(-3.0905753) q[0];
sx q[0];
rz(-1.464123) q[0];
rz(-1.4785712) q[1];
sx q[1];
rz(-1.1590978) q[1];
sx q[1];
rz(1.2082072) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0776652) q[0];
sx q[0];
rz(-2.0443235) q[0];
sx q[0];
rz(-2.5350476) q[0];
x q[1];
rz(1.1651498) q[2];
sx q[2];
rz(-2.6313836) q[2];
sx q[2];
rz(-2.0070397) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3390926) q[1];
sx q[1];
rz(-2.8475757) q[1];
sx q[1];
rz(1.1400998) q[1];
rz(-pi) q[2];
rz(-0.34388108) q[3];
sx q[3];
rz(-2.6626427) q[3];
sx q[3];
rz(-2.6987181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6055484) q[2];
sx q[2];
rz(-0.41993419) q[2];
sx q[2];
rz(1.4844683) q[2];
rz(-3.1260955) q[3];
sx q[3];
rz(-1.213538) q[3];
sx q[3];
rz(1.9626455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.991796) q[0];
sx q[0];
rz(-1.3301671) q[0];
sx q[0];
rz(0.27161828) q[0];
rz(2.244921) q[1];
sx q[1];
rz(-2.6397557) q[1];
sx q[1];
rz(1.2976049) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5974242) q[0];
sx q[0];
rz(-2.6784458) q[0];
sx q[0];
rz(-1.5933655) q[0];
x q[1];
rz(-2.5141352) q[2];
sx q[2];
rz(-0.59840032) q[2];
sx q[2];
rz(2.5638169) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1101164) q[1];
sx q[1];
rz(-1.3240485) q[1];
sx q[1];
rz(-0.11746789) q[1];
rz(-pi) q[2];
rz(-0.0053169189) q[3];
sx q[3];
rz(-2.379619) q[3];
sx q[3];
rz(-3.0865106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4012332) q[2];
sx q[2];
rz(-2.3637502) q[2];
sx q[2];
rz(-0.05376251) q[2];
rz(1.8476123) q[3];
sx q[3];
rz(-0.79837489) q[3];
sx q[3];
rz(0.23538858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2006705) q[0];
sx q[0];
rz(-0.5680474) q[0];
sx q[0];
rz(-1.6773552) q[0];
rz(-2.0109406) q[1];
sx q[1];
rz(-2.407275) q[1];
sx q[1];
rz(0.11437036) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8371724) q[0];
sx q[0];
rz(-0.74203287) q[0];
sx q[0];
rz(-2.4993308) q[0];
x q[1];
rz(1.0724154) q[2];
sx q[2];
rz(-0.97626462) q[2];
sx q[2];
rz(1.6785113) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3721766) q[1];
sx q[1];
rz(-1.5443364) q[1];
sx q[1];
rz(-0.89511223) q[1];
rz(-1.5773729) q[3];
sx q[3];
rz(-1.3727756) q[3];
sx q[3];
rz(1.1634367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4734681) q[2];
sx q[2];
rz(-1.5359842) q[2];
sx q[2];
rz(-3.065897) q[2];
rz(-2.7068052) q[3];
sx q[3];
rz(-1.9861168) q[3];
sx q[3];
rz(-3.0619612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39744034) q[0];
sx q[0];
rz(-1.8295153) q[0];
sx q[0];
rz(-2.721526) q[0];
rz(-2.9282667) q[1];
sx q[1];
rz(-1.408564) q[1];
sx q[1];
rz(1.2423645) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38753375) q[0];
sx q[0];
rz(-2.2926169) q[0];
sx q[0];
rz(2.5639064) q[0];
rz(1.2790643) q[2];
sx q[2];
rz(-1.1559249) q[2];
sx q[2];
rz(1.3049098) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1429726) q[1];
sx q[1];
rz(-3.0398453) q[1];
sx q[1];
rz(1.7228026) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.20988864) q[3];
sx q[3];
rz(-2.915463) q[3];
sx q[3];
rz(-1.9551639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3199978) q[2];
sx q[2];
rz(-2.2264806) q[2];
sx q[2];
rz(-0.37459174) q[2];
rz(-2.1290667) q[3];
sx q[3];
rz(-0.39207021) q[3];
sx q[3];
rz(-2.4969126) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11151611) q[0];
sx q[0];
rz(-2.6426297) q[0];
sx q[0];
rz(2.7110355) q[0];
rz(2.997609) q[1];
sx q[1];
rz(-1.0824243) q[1];
sx q[1];
rz(2.5103501) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3595561) q[0];
sx q[0];
rz(-2.5304473) q[0];
sx q[0];
rz(3.0474328) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9955533) q[2];
sx q[2];
rz(-1.1542873) q[2];
sx q[2];
rz(-3.0768968) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5938607) q[1];
sx q[1];
rz(-2.0640848) q[1];
sx q[1];
rz(0.6823205) q[1];
rz(-2.1926375) q[3];
sx q[3];
rz(-2.6027711) q[3];
sx q[3];
rz(-2.1486189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.97079078) q[2];
sx q[2];
rz(-2.0826714) q[2];
sx q[2];
rz(-1.4405174) q[2];
rz(1.3614281) q[3];
sx q[3];
rz(-1.1037339) q[3];
sx q[3];
rz(0.48719278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.454527) q[0];
sx q[0];
rz(-0.28959689) q[0];
sx q[0];
rz(0.92426306) q[0];
rz(-2.848792) q[1];
sx q[1];
rz(-1.2288789) q[1];
sx q[1];
rz(-0.80088314) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7445114) q[0];
sx q[0];
rz(-1.9430706) q[0];
sx q[0];
rz(2.8651994) q[0];
rz(-2.7960294) q[2];
sx q[2];
rz(-0.93831944) q[2];
sx q[2];
rz(-1.3407624) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.85383254) q[1];
sx q[1];
rz(-1.6759071) q[1];
sx q[1];
rz(-1.5413324) q[1];
x q[2];
rz(3.1195779) q[3];
sx q[3];
rz(-1.5436633) q[3];
sx q[3];
rz(2.5282945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.822829) q[2];
sx q[2];
rz(-1.2206565) q[2];
sx q[2];
rz(1.8611543) q[2];
rz(-0.66796962) q[3];
sx q[3];
rz(-0.48798713) q[3];
sx q[3];
rz(0.41332301) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8845344) q[0];
sx q[0];
rz(-1.9030544) q[0];
sx q[0];
rz(1.1827693) q[0];
rz(0.23712748) q[1];
sx q[1];
rz(-3.0393937) q[1];
sx q[1];
rz(0.059344083) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5049669) q[0];
sx q[0];
rz(-1.6904181) q[0];
sx q[0];
rz(-1.8455532) q[0];
x q[1];
rz(-0.72513942) q[2];
sx q[2];
rz(-1.983101) q[2];
sx q[2];
rz(-0.39620846) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0036111) q[1];
sx q[1];
rz(-1.4661769) q[1];
sx q[1];
rz(1.8358747) q[1];
x q[2];
rz(-3.0640934) q[3];
sx q[3];
rz(-1.2480253) q[3];
sx q[3];
rz(-0.68068824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5817029) q[2];
sx q[2];
rz(-0.54055944) q[2];
sx q[2];
rz(1.7891368) q[2];
rz(-1.3672359) q[3];
sx q[3];
rz(-1.4227899) q[3];
sx q[3];
rz(-2.2860315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2818114) q[0];
sx q[0];
rz(-0.78524041) q[0];
sx q[0];
rz(-0.45968858) q[0];
rz(-3.1135318) q[1];
sx q[1];
rz(-1.1496081) q[1];
sx q[1];
rz(1.185816) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5494649) q[0];
sx q[0];
rz(-2.87232) q[0];
sx q[0];
rz(0.16945355) q[0];
rz(1.6721647) q[2];
sx q[2];
rz(-1.4918054) q[2];
sx q[2];
rz(1.7192507) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8227905) q[1];
sx q[1];
rz(-1.0981961) q[1];
sx q[1];
rz(2.5986555) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7544153) q[3];
sx q[3];
rz(-2.4868591) q[3];
sx q[3];
rz(2.2019405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6431553) q[2];
sx q[2];
rz(-2.2648621) q[2];
sx q[2];
rz(-0.15677491) q[2];
rz(2.5458941) q[3];
sx q[3];
rz(-2.7955293) q[3];
sx q[3];
rz(-3.116385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51782411) q[0];
sx q[0];
rz(-2.1110004) q[0];
sx q[0];
rz(0.18381707) q[0];
rz(3.0635762) q[1];
sx q[1];
rz(-1.0573496) q[1];
sx q[1];
rz(-2.877291) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.030609) q[0];
sx q[0];
rz(-1.8389337) q[0];
sx q[0];
rz(-0.58695745) q[0];
rz(-pi) q[1];
x q[1];
rz(1.043522) q[2];
sx q[2];
rz(-2.1855178) q[2];
sx q[2];
rz(-1.5103024) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.42929998) q[1];
sx q[1];
rz(-1.7279062) q[1];
sx q[1];
rz(2.9595441) q[1];
rz(-2.9654042) q[3];
sx q[3];
rz(-2.6399351) q[3];
sx q[3];
rz(2.9210726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.8514303) q[2];
sx q[2];
rz(-2.1992079) q[2];
sx q[2];
rz(-2.9128722) q[2];
rz(-1.4043572) q[3];
sx q[3];
rz(-2.2018933) q[3];
sx q[3];
rz(-3.0586045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1514773) q[0];
sx q[0];
rz(-2.293312) q[0];
sx q[0];
rz(1.5207416) q[0];
rz(-0.71113853) q[1];
sx q[1];
rz(-1.3093206) q[1];
sx q[1];
rz(0.73285229) q[1];
rz(-1.0719094) q[2];
sx q[2];
rz(-1.5317393) q[2];
sx q[2];
rz(-2.692937) q[2];
rz(-1.0300954) q[3];
sx q[3];
rz(-1.0813011) q[3];
sx q[3];
rz(-2.4745221) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
