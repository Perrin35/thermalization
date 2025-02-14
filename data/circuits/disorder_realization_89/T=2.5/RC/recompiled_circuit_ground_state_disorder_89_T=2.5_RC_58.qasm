OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9417579) q[0];
sx q[0];
rz(6.9747337) q[0];
sx q[0];
rz(11.792948) q[0];
rz(3.050488) q[1];
sx q[1];
rz(-0.79371047) q[1];
sx q[1];
rz(-1.33338) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7156336) q[0];
sx q[0];
rz(-1.5797401) q[0];
sx q[0];
rz(-2.483213) q[0];
rz(-1.0716419) q[2];
sx q[2];
rz(-1.7737651) q[2];
sx q[2];
rz(-2.2405193) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4845407) q[1];
sx q[1];
rz(-2.1738677) q[1];
sx q[1];
rz(-2.6588427) q[1];
rz(-1.638014) q[3];
sx q[3];
rz(-1.5537694) q[3];
sx q[3];
rz(-2.5650781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.29204631) q[2];
sx q[2];
rz(-1.8201733) q[2];
sx q[2];
rz(-2.7267961) q[2];
rz(-1.8997806) q[3];
sx q[3];
rz(-1.1153699) q[3];
sx q[3];
rz(2.950086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6280129) q[0];
sx q[0];
rz(-0.11592557) q[0];
sx q[0];
rz(-0.43981788) q[0];
rz(-0.10305931) q[1];
sx q[1];
rz(-0.28918806) q[1];
sx q[1];
rz(1.1691079) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3713908) q[0];
sx q[0];
rz(-1.4318716) q[0];
sx q[0];
rz(2.2004805) q[0];
x q[1];
rz(-0.051889433) q[2];
sx q[2];
rz(-0.98116335) q[2];
sx q[2];
rz(2.6419214) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4354187) q[1];
sx q[1];
rz(-1.7507554) q[1];
sx q[1];
rz(0.91689308) q[1];
rz(1.7984932) q[3];
sx q[3];
rz(-2.4191354) q[3];
sx q[3];
rz(-1.8060773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.49587265) q[2];
sx q[2];
rz(-1.7503259) q[2];
sx q[2];
rz(0.24620852) q[2];
rz(-1.5919033) q[3];
sx q[3];
rz(-2.527745) q[3];
sx q[3];
rz(-1.7483819) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.15568) q[0];
sx q[0];
rz(-1.14013) q[0];
sx q[0];
rz(0.21173665) q[0];
rz(0.011292975) q[1];
sx q[1];
rz(-0.79835049) q[1];
sx q[1];
rz(1.4439772) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2560476) q[0];
sx q[0];
rz(-3.1064337) q[0];
sx q[0];
rz(3.1395182) q[0];
x q[1];
rz(-1.774774) q[2];
sx q[2];
rz(-0.24572554) q[2];
sx q[2];
rz(-1.7468921) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9992175) q[1];
sx q[1];
rz(-2.1793167) q[1];
sx q[1];
rz(-1.1097679) q[1];
rz(-pi) q[2];
rz(2.9559198) q[3];
sx q[3];
rz(-0.45760307) q[3];
sx q[3];
rz(2.5724907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.4766562) q[2];
sx q[2];
rz(-2.0659451) q[2];
sx q[2];
rz(-0.54534379) q[2];
rz(2.2319345) q[3];
sx q[3];
rz(-2.6792512) q[3];
sx q[3];
rz(-1.9285412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42105168) q[0];
sx q[0];
rz(-2.4626829) q[0];
sx q[0];
rz(-2.3140267) q[0];
rz(1.958581) q[1];
sx q[1];
rz(-1.8667826) q[1];
sx q[1];
rz(1.8756728) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8592081) q[0];
sx q[0];
rz(-3.0050238) q[0];
sx q[0];
rz(2.300416) q[0];
rz(-pi) q[1];
rz(-1.5180196) q[2];
sx q[2];
rz(-0.20773331) q[2];
sx q[2];
rz(-0.29333255) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3490679) q[1];
sx q[1];
rz(-1.6025839) q[1];
sx q[1];
rz(-2.1963901) q[1];
rz(1.3807348) q[3];
sx q[3];
rz(-1.7281088) q[3];
sx q[3];
rz(-0.33657956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8852692) q[2];
sx q[2];
rz(-1.4726535) q[2];
sx q[2];
rz(3.0598158) q[2];
rz(0.30803382) q[3];
sx q[3];
rz(-0.48397288) q[3];
sx q[3];
rz(-1.5380194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1643739) q[0];
sx q[0];
rz(-2.866221) q[0];
sx q[0];
rz(-0.070505738) q[0];
rz(-0.33440691) q[1];
sx q[1];
rz(-2.6102378) q[1];
sx q[1];
rz(-2.6883584) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95246661) q[0];
sx q[0];
rz(-2.1299908) q[0];
sx q[0];
rz(-2.1711012) q[0];
rz(-pi) q[1];
rz(1.9627175) q[2];
sx q[2];
rz(-1.3597314) q[2];
sx q[2];
rz(-3.0818194) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0754521) q[1];
sx q[1];
rz(-1.4601216) q[1];
sx q[1];
rz(2.8933725) q[1];
rz(2.2599561) q[3];
sx q[3];
rz(-0.83587468) q[3];
sx q[3];
rz(0.78208941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3004904) q[2];
sx q[2];
rz(-1.6359685) q[2];
sx q[2];
rz(2.9312768) q[2];
rz(-0.39120832) q[3];
sx q[3];
rz(-2.936383) q[3];
sx q[3];
rz(2.6989663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6415569) q[0];
sx q[0];
rz(-2.0099202) q[0];
sx q[0];
rz(-2.9697707) q[0];
rz(-1.3459282) q[1];
sx q[1];
rz(-2.185952) q[1];
sx q[1];
rz(-1.1489493) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1412107) q[0];
sx q[0];
rz(-0.11185574) q[0];
sx q[0];
rz(-3.0684708) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23229084) q[2];
sx q[2];
rz(-1.8477167) q[2];
sx q[2];
rz(-2.0611219) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2021804) q[1];
sx q[1];
rz(-1.3883852) q[1];
sx q[1];
rz(-0.98947452) q[1];
rz(-pi) q[2];
rz(1.9257718) q[3];
sx q[3];
rz(-2.6067511) q[3];
sx q[3];
rz(-2.8876784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.12080869) q[2];
sx q[2];
rz(-0.9641996) q[2];
sx q[2];
rz(2.7210893) q[2];
rz(1.486982) q[3];
sx q[3];
rz(-1.9872811) q[3];
sx q[3];
rz(1.1544863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7789975) q[0];
sx q[0];
rz(-1.2263466) q[0];
sx q[0];
rz(2.763789) q[0];
rz(-1.3899577) q[1];
sx q[1];
rz(-1.4327587) q[1];
sx q[1];
rz(-0.43209824) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88277757) q[0];
sx q[0];
rz(-1.7452137) q[0];
sx q[0];
rz(2.9091878) q[0];
rz(2.7977706) q[2];
sx q[2];
rz(-3.0431586) q[2];
sx q[2];
rz(-1.4296895) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3258265) q[1];
sx q[1];
rz(-1.6066505) q[1];
sx q[1];
rz(0.46200606) q[1];
rz(-2.7338846) q[3];
sx q[3];
rz(-2.2033436) q[3];
sx q[3];
rz(-0.11128259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3720588) q[2];
sx q[2];
rz(-0.62166119) q[2];
sx q[2];
rz(0.2335693) q[2];
rz(1.6893049) q[3];
sx q[3];
rz(-1.7100311) q[3];
sx q[3];
rz(1.5610032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5772783) q[0];
sx q[0];
rz(-1.0294585) q[0];
sx q[0];
rz(0.3197864) q[0];
rz(-0.86209595) q[1];
sx q[1];
rz(-1.2513221) q[1];
sx q[1];
rz(1.7822942) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0523543) q[0];
sx q[0];
rz(-2.1349094) q[0];
sx q[0];
rz(-0.67676877) q[0];
x q[1];
rz(0.14131693) q[2];
sx q[2];
rz(-1.769678) q[2];
sx q[2];
rz(-0.74706739) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.055147) q[1];
sx q[1];
rz(-2.3999955) q[1];
sx q[1];
rz(-2.8838709) q[1];
rz(-pi) q[2];
rz(-1.9407746) q[3];
sx q[3];
rz(-0.66744679) q[3];
sx q[3];
rz(-0.98946179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.84997815) q[2];
sx q[2];
rz(-2.7225967) q[2];
sx q[2];
rz(-0.46763793) q[2];
rz(2.6575798) q[3];
sx q[3];
rz(-1.368112) q[3];
sx q[3];
rz(1.8638994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(-0.17396946) q[0];
sx q[0];
rz(-1.4894217) q[0];
sx q[0];
rz(1.1349348) q[0];
rz(3.0813772) q[1];
sx q[1];
rz(-0.73679149) q[1];
sx q[1];
rz(1.6103475) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9519015) q[0];
sx q[0];
rz(-0.78329059) q[0];
sx q[0];
rz(-0.89166151) q[0];
rz(-1.0410454) q[2];
sx q[2];
rz(-2.1868631) q[2];
sx q[2];
rz(3.0977647) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.169226) q[1];
sx q[1];
rz(-1.3080096) q[1];
sx q[1];
rz(-1.4701094) q[1];
rz(0.78218145) q[3];
sx q[3];
rz(-1.5429268) q[3];
sx q[3];
rz(-3.0620344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.81426364) q[2];
sx q[2];
rz(-1.9244104) q[2];
sx q[2];
rz(0.72861707) q[2];
rz(-0.21555756) q[3];
sx q[3];
rz(-1.7451127) q[3];
sx q[3];
rz(2.047915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2586867) q[0];
sx q[0];
rz(-0.031143324) q[0];
sx q[0];
rz(2.8236142) q[0];
rz(1.7440354) q[1];
sx q[1];
rz(-1.8122858) q[1];
sx q[1];
rz(2.8584282) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69071373) q[0];
sx q[0];
rz(-2.9895999) q[0];
sx q[0];
rz(2.1426523) q[0];
x q[1];
rz(0.0816143) q[2];
sx q[2];
rz(-2.3858262) q[2];
sx q[2];
rz(2.4163742) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4589112) q[1];
sx q[1];
rz(-1.1106756) q[1];
sx q[1];
rz(0.7724055) q[1];
rz(1.3628325) q[3];
sx q[3];
rz(-1.4803518) q[3];
sx q[3];
rz(-1.2048723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2240923) q[2];
sx q[2];
rz(-1.516569) q[2];
sx q[2];
rz(3.0523114) q[2];
rz(-1.8381522) q[3];
sx q[3];
rz(-0.83804122) q[3];
sx q[3];
rz(-2.1681521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3598809) q[0];
sx q[0];
rz(-0.64890535) q[0];
sx q[0];
rz(-2.1606408) q[0];
rz(2.5557062) q[1];
sx q[1];
rz(-1.8460907) q[1];
sx q[1];
rz(1.5150217) q[1];
rz(2.9286251) q[2];
sx q[2];
rz(-1.3942547) q[2];
sx q[2];
rz(-0.19993776) q[2];
rz(-1.951915) q[3];
sx q[3];
rz(-1.3314691) q[3];
sx q[3];
rz(-0.21503147) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
