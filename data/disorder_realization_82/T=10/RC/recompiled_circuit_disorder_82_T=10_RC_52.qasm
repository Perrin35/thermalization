OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80469552) q[0];
sx q[0];
rz(-1.0372294) q[0];
sx q[0];
rz(-2.7859935) q[0];
rz(2.8388677) q[1];
sx q[1];
rz(-1.0441138) q[1];
sx q[1];
rz(1.27966) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.049392603) q[0];
sx q[0];
rz(-2.85673) q[0];
sx q[0];
rz(-1.0546791) q[0];
rz(-1.3470634) q[2];
sx q[2];
rz(-1.9828084) q[2];
sx q[2];
rz(1.7681233) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0496088) q[1];
sx q[1];
rz(-2.1503452) q[1];
sx q[1];
rz(-2.0642573) q[1];
x q[2];
rz(1.3650595) q[3];
sx q[3];
rz(-1.7997026) q[3];
sx q[3];
rz(-2.8347901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1352284) q[2];
sx q[2];
rz(-2.2029115) q[2];
sx q[2];
rz(0.65650666) q[2];
rz(-2.4025829) q[3];
sx q[3];
rz(-0.46027547) q[3];
sx q[3];
rz(0.41729331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0083369) q[0];
sx q[0];
rz(-2.3454741) q[0];
sx q[0];
rz(0.59536368) q[0];
rz(3.0796675) q[1];
sx q[1];
rz(-1.9134816) q[1];
sx q[1];
rz(2.6541236) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9092642) q[0];
sx q[0];
rz(-3.0380913) q[0];
sx q[0];
rz(-0.89719015) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42178085) q[2];
sx q[2];
rz(-1.9739208) q[2];
sx q[2];
rz(0.031907206) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0812379) q[1];
sx q[1];
rz(-1.758467) q[1];
sx q[1];
rz(0.16138046) q[1];
rz(-pi) q[2];
rz(-3.0549166) q[3];
sx q[3];
rz(-1.254734) q[3];
sx q[3];
rz(2.812127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9618824) q[2];
sx q[2];
rz(-1.8790745) q[2];
sx q[2];
rz(2.7462192) q[2];
rz(-2.0987434) q[3];
sx q[3];
rz(-2.6104749) q[3];
sx q[3];
rz(0.038671967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19668002) q[0];
sx q[0];
rz(-2.0356464) q[0];
sx q[0];
rz(-0.19038598) q[0];
rz(0.12292513) q[1];
sx q[1];
rz(-0.38750896) q[1];
sx q[1];
rz(0.22274676) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4653141) q[0];
sx q[0];
rz(-1.9440117) q[0];
sx q[0];
rz(-0.42155427) q[0];
rz(0.22521714) q[2];
sx q[2];
rz(-1.2080964) q[2];
sx q[2];
rz(-0.97937102) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60274117) q[1];
sx q[1];
rz(-0.95898421) q[1];
sx q[1];
rz(-1.6536504) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1915216) q[3];
sx q[3];
rz(-0.40682236) q[3];
sx q[3];
rz(1.0328968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2591851) q[2];
sx q[2];
rz(-1.8683445) q[2];
sx q[2];
rz(1.1941236) q[2];
rz(-1.0549226) q[3];
sx q[3];
rz(-1.4849562) q[3];
sx q[3];
rz(0.96364337) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7524183) q[0];
sx q[0];
rz(-0.27802935) q[0];
sx q[0];
rz(1.4032723) q[0];
rz(-1.4933043) q[1];
sx q[1];
rz(-1.5999258) q[1];
sx q[1];
rz(2.2213675) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7677778) q[0];
sx q[0];
rz(-2.2245363) q[0];
sx q[0];
rz(-0.52315229) q[0];
x q[1];
rz(0.36107365) q[2];
sx q[2];
rz(-1.4898584) q[2];
sx q[2];
rz(-2.0780448) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8060382) q[1];
sx q[1];
rz(-1.1127377) q[1];
sx q[1];
rz(2.4469417) q[1];
rz(1.7761049) q[3];
sx q[3];
rz(-1.218154) q[3];
sx q[3];
rz(2.7384788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9437287) q[2];
sx q[2];
rz(-1.7284164) q[2];
sx q[2];
rz(-1.2801923) q[2];
rz(2.3222893) q[3];
sx q[3];
rz(-1.8236482) q[3];
sx q[3];
rz(1.7839446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.458805) q[0];
sx q[0];
rz(-1.7651432) q[0];
sx q[0];
rz(1.4047594) q[0];
rz(0.76830307) q[1];
sx q[1];
rz(-0.65015018) q[1];
sx q[1];
rz(-0.4531025) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83349193) q[0];
sx q[0];
rz(-1.6080603) q[0];
sx q[0];
rz(-1.4971855) q[0];
rz(-pi) q[1];
x q[1];
rz(0.034899072) q[2];
sx q[2];
rz(-1.7136095) q[2];
sx q[2];
rz(2.3226483) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0616152) q[1];
sx q[1];
rz(-1.4887759) q[1];
sx q[1];
rz(1.7916937) q[1];
rz(-pi) q[2];
rz(2.1702607) q[3];
sx q[3];
rz(-1.9595651) q[3];
sx q[3];
rz(-0.99861162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.12864628) q[2];
sx q[2];
rz(-0.83798989) q[2];
sx q[2];
rz(-1.6298693) q[2];
rz(2.417918) q[3];
sx q[3];
rz(-1.2440163) q[3];
sx q[3];
rz(-0.85030142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033427514) q[0];
sx q[0];
rz(-1.4135679) q[0];
sx q[0];
rz(2.9034555) q[0];
rz(-0.37995964) q[1];
sx q[1];
rz(-1.0521051) q[1];
sx q[1];
rz(-1.628081) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1171922) q[0];
sx q[0];
rz(-1.8937366) q[0];
sx q[0];
rz(0.49844235) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0797708) q[2];
sx q[2];
rz(-1.1218058) q[2];
sx q[2];
rz(-2.2668554) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7442419) q[1];
sx q[1];
rz(-0.89351082) q[1];
sx q[1];
rz(0.97126295) q[1];
rz(-pi) q[2];
rz(2.7558277) q[3];
sx q[3];
rz(-1.9915951) q[3];
sx q[3];
rz(-2.5667131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0825519) q[2];
sx q[2];
rz(-2.5566792) q[2];
sx q[2];
rz(1.1435821) q[2];
rz(3.011076) q[3];
sx q[3];
rz(-1.427622) q[3];
sx q[3];
rz(-0.17091621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1259574) q[0];
sx q[0];
rz(-2.1645249) q[0];
sx q[0];
rz(-2.7440199) q[0];
rz(1.827084) q[1];
sx q[1];
rz(-0.54388261) q[1];
sx q[1];
rz(-1.1869173) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8966973) q[0];
sx q[0];
rz(-0.58422409) q[0];
sx q[0];
rz(-1.0446192) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45849623) q[2];
sx q[2];
rz(-0.065932238) q[2];
sx q[2];
rz(2.3022848) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.742332) q[1];
sx q[1];
rz(-0.84335828) q[1];
sx q[1];
rz(2.9620693) q[1];
rz(-2.5116634) q[3];
sx q[3];
rz(-1.3327206) q[3];
sx q[3];
rz(1.5671974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.879803) q[2];
sx q[2];
rz(-1.3061085) q[2];
sx q[2];
rz(1.3558033) q[2];
rz(1.5073744) q[3];
sx q[3];
rz(-2.5528788) q[3];
sx q[3];
rz(0.9986977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3223406) q[0];
sx q[0];
rz(-1.9544019) q[0];
sx q[0];
rz(0.85574714) q[0];
rz(0.021727173) q[1];
sx q[1];
rz(-2.1179492) q[1];
sx q[1];
rz(-1.0303248) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.94901) q[0];
sx q[0];
rz(-2.2826676) q[0];
sx q[0];
rz(2.683995) q[0];
x q[1];
rz(-1.0439992) q[2];
sx q[2];
rz(-1.6128522) q[2];
sx q[2];
rz(-0.85925697) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1134539) q[1];
sx q[1];
rz(-0.33543643) q[1];
sx q[1];
rz(-1.0949275) q[1];
rz(-pi) q[2];
rz(0.52328531) q[3];
sx q[3];
rz(-2.0998294) q[3];
sx q[3];
rz(0.39213359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8490303) q[2];
sx q[2];
rz(-1.2529255) q[2];
sx q[2];
rz(-0.070177468) q[2];
rz(-0.82693806) q[3];
sx q[3];
rz(-1.6707872) q[3];
sx q[3];
rz(0.58644811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4940015) q[0];
sx q[0];
rz(-1.4275455) q[0];
sx q[0];
rz(0.31627396) q[0];
rz(1.0519741) q[1];
sx q[1];
rz(-0.72967044) q[1];
sx q[1];
rz(2.0064328) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90056706) q[0];
sx q[0];
rz(-1.3824438) q[0];
sx q[0];
rz(-2.7450949) q[0];
rz(-pi) q[1];
rz(-1.4831545) q[2];
sx q[2];
rz(-0.32155514) q[2];
sx q[2];
rz(1.1631539) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4371722) q[1];
sx q[1];
rz(-1.6295027) q[1];
sx q[1];
rz(-1.0874332) q[1];
rz(-pi) q[2];
rz(1.5376066) q[3];
sx q[3];
rz(-0.6233218) q[3];
sx q[3];
rz(-0.25191307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6241374) q[2];
sx q[2];
rz(-0.74666658) q[2];
sx q[2];
rz(-1.5448145) q[2];
rz(2.4638713) q[3];
sx q[3];
rz(-2.2993408) q[3];
sx q[3];
rz(-2.8519582) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.233376) q[0];
sx q[0];
rz(-2.4514618) q[0];
sx q[0];
rz(-0.35183516) q[0];
rz(2.8219163) q[1];
sx q[1];
rz(-0.37477481) q[1];
sx q[1];
rz(2.9454254) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2575532) q[0];
sx q[0];
rz(-1.0172052) q[0];
sx q[0];
rz(-2.2474225) q[0];
rz(1.9253795) q[2];
sx q[2];
rz(-2.9973929) q[2];
sx q[2];
rz(-1.2186288) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0077121) q[1];
sx q[1];
rz(-2.4073283) q[1];
sx q[1];
rz(-1.5937362) q[1];
x q[2];
rz(2.4862643) q[3];
sx q[3];
rz(-2.1851551) q[3];
sx q[3];
rz(-0.67656803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3907884) q[2];
sx q[2];
rz(-1.6395586) q[2];
sx q[2];
rz(-3.0604559) q[2];
rz(1.9598512) q[3];
sx q[3];
rz(-2.2406082) q[3];
sx q[3];
rz(-1.0749764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1186196) q[0];
sx q[0];
rz(-1.5587627) q[0];
sx q[0];
rz(-1.8297304) q[0];
rz(-1.3148057) q[1];
sx q[1];
rz(-0.79827764) q[1];
sx q[1];
rz(-1.6006443) q[1];
rz(2.9076004) q[2];
sx q[2];
rz(-1.5237612) q[2];
sx q[2];
rz(0.16618726) q[2];
rz(-1.374791) q[3];
sx q[3];
rz(-0.23402611) q[3];
sx q[3];
rz(0.2454161) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
