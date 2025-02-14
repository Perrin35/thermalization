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
rz(-1.8974798) q[0];
sx q[0];
rz(-0.31960684) q[0];
sx q[0];
rz(-0.50080103) q[0];
rz(-2.8473941) q[1];
sx q[1];
rz(-0.74217141) q[1];
sx q[1];
rz(2.5941526) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8590605) q[0];
sx q[0];
rz(-1.7164537) q[0];
sx q[0];
rz(0.3126302) q[0];
rz(-pi) q[1];
rz(0.51551657) q[2];
sx q[2];
rz(-0.52243865) q[2];
sx q[2];
rz(-1.9439657) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.77078041) q[1];
sx q[1];
rz(-0.51518422) q[1];
sx q[1];
rz(2.9574802) q[1];
rz(-2.7846401) q[3];
sx q[3];
rz(-1.0610749) q[3];
sx q[3];
rz(-1.3689643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1529634) q[2];
sx q[2];
rz(-2.2965778) q[2];
sx q[2];
rz(0.60626924) q[2];
rz(1.2074977) q[3];
sx q[3];
rz(-1.2946125) q[3];
sx q[3];
rz(-1.5906364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2400874) q[0];
sx q[0];
rz(-1.5272239) q[0];
sx q[0];
rz(-0.69183451) q[0];
rz(1.8531307) q[1];
sx q[1];
rz(-1.9970147) q[1];
sx q[1];
rz(2.8537234) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26507227) q[0];
sx q[0];
rz(-1.2848043) q[0];
sx q[0];
rz(0.13354451) q[0];
x q[1];
rz(-0.68142724) q[2];
sx q[2];
rz(-2.3286901) q[2];
sx q[2];
rz(-2.5035448) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.58093666) q[1];
sx q[1];
rz(-1.1252778) q[1];
sx q[1];
rz(-1.5195283) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4631264) q[3];
sx q[3];
rz(-1.7070908) q[3];
sx q[3];
rz(-1.8267269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79491478) q[2];
sx q[2];
rz(-1.4823464) q[2];
sx q[2];
rz(2.1686926) q[2];
rz(1.8074624) q[3];
sx q[3];
rz(-2.3502246) q[3];
sx q[3];
rz(2.5202675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36488229) q[0];
sx q[0];
rz(-2.8948687) q[0];
sx q[0];
rz(2.8926988) q[0];
rz(-1.4502672) q[1];
sx q[1];
rz(-1.9906809) q[1];
sx q[1];
rz(0.90363622) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2657651) q[0];
sx q[0];
rz(-0.96185124) q[0];
sx q[0];
rz(1.5399571) q[0];
x q[1];
rz(1.2799706) q[2];
sx q[2];
rz(-2.7827459) q[2];
sx q[2];
rz(-1.0076696) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.98164046) q[1];
sx q[1];
rz(-2.2800696) q[1];
sx q[1];
rz(-1.0538799) q[1];
rz(0.60402212) q[3];
sx q[3];
rz(-1.9088749) q[3];
sx q[3];
rz(0.018751831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9768208) q[2];
sx q[2];
rz(-2.4103319) q[2];
sx q[2];
rz(-2.8110647) q[2];
rz(2.9183689) q[3];
sx q[3];
rz(-1.1130755) q[3];
sx q[3];
rz(2.7631675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44499236) q[0];
sx q[0];
rz(-1.7855423) q[0];
sx q[0];
rz(-2.9997605) q[0];
rz(0.089135535) q[1];
sx q[1];
rz(-2.8991957) q[1];
sx q[1];
rz(0.95318046) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0713816) q[0];
sx q[0];
rz(-1.9527023) q[0];
sx q[0];
rz(-0.75526636) q[0];
rz(-pi) q[1];
rz(2.761327) q[2];
sx q[2];
rz(-1.6820246) q[2];
sx q[2];
rz(-0.56651543) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.64884463) q[1];
sx q[1];
rz(-2.358583) q[1];
sx q[1];
rz(-2.3451508) q[1];
rz(0.35928407) q[3];
sx q[3];
rz(-2.6169638) q[3];
sx q[3];
rz(-2.5398162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3700833) q[2];
sx q[2];
rz(-0.13088317) q[2];
sx q[2];
rz(0.5292325) q[2];
rz(-1.0445163) q[3];
sx q[3];
rz(-1.3890725) q[3];
sx q[3];
rz(0.15531003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20045497) q[0];
sx q[0];
rz(-1.8998242) q[0];
sx q[0];
rz(-0.45503765) q[0];
rz(3.1269238) q[1];
sx q[1];
rz(-0.55611062) q[1];
sx q[1];
rz(-0.98247772) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9519314) q[0];
sx q[0];
rz(-1.8704544) q[0];
sx q[0];
rz(1.4656187) q[0];
rz(-pi) q[1];
rz(0.90679866) q[2];
sx q[2];
rz(-2.325146) q[2];
sx q[2];
rz(1.0524257) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8019218) q[1];
sx q[1];
rz(-1.675444) q[1];
sx q[1];
rz(0.27709477) q[1];
x q[2];
rz(-1.807688) q[3];
sx q[3];
rz(-1.8474527) q[3];
sx q[3];
rz(1.93678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.602953) q[2];
sx q[2];
rz(-2.2621138) q[2];
sx q[2];
rz(-2.5243916) q[2];
rz(-1.7547102) q[3];
sx q[3];
rz(-2.230547) q[3];
sx q[3];
rz(2.9673747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0844326) q[0];
sx q[0];
rz(-0.76157695) q[0];
sx q[0];
rz(2.1630951) q[0];
rz(1.017336) q[1];
sx q[1];
rz(-2.8637736) q[1];
sx q[1];
rz(0.4496347) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.250424) q[0];
sx q[0];
rz(-0.40968597) q[0];
sx q[0];
rz(2.4158359) q[0];
x q[1];
rz(-1.1005747) q[2];
sx q[2];
rz(-0.82056773) q[2];
sx q[2];
rz(2.5448397) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2404473) q[1];
sx q[1];
rz(-2.6679651) q[1];
sx q[1];
rz(1.4722826) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1862667) q[3];
sx q[3];
rz(-0.98444533) q[3];
sx q[3];
rz(2.9406656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6189239) q[2];
sx q[2];
rz(-1.7974682) q[2];
sx q[2];
rz(1.3783003) q[2];
rz(-1.1528692) q[3];
sx q[3];
rz(-1.1533573) q[3];
sx q[3];
rz(2.8809179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5744837) q[0];
sx q[0];
rz(-2.4686047) q[0];
sx q[0];
rz(-2.8048977) q[0];
rz(-0.7252655) q[1];
sx q[1];
rz(-1.5545574) q[1];
sx q[1];
rz(-0.38549647) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9850901) q[0];
sx q[0];
rz(-0.96961951) q[0];
sx q[0];
rz(-1.9999835) q[0];
x q[1];
rz(-1.8202838) q[2];
sx q[2];
rz(-1.7091342) q[2];
sx q[2];
rz(2.4398566) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.817877) q[1];
sx q[1];
rz(-2.1855832) q[1];
sx q[1];
rz(-0.35732689) q[1];
rz(1.9970421) q[3];
sx q[3];
rz(-0.9916456) q[3];
sx q[3];
rz(0.338011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5943299) q[2];
sx q[2];
rz(-1.9219834) q[2];
sx q[2];
rz(-2.8426113) q[2];
rz(2.3494521) q[3];
sx q[3];
rz(-2.7463089) q[3];
sx q[3];
rz(1.7637174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2250273) q[0];
sx q[0];
rz(-1.7074317) q[0];
sx q[0];
rz(0.7780956) q[0];
rz(-1.3968702) q[1];
sx q[1];
rz(-0.94051802) q[1];
sx q[1];
rz(-3.0527589) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7185811) q[0];
sx q[0];
rz(-1.4664343) q[0];
sx q[0];
rz(0.34684885) q[0];
rz(-pi) q[1];
rz(-0.93508945) q[2];
sx q[2];
rz(-2.3176821) q[2];
sx q[2];
rz(-1.7721792) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.962304) q[1];
sx q[1];
rz(-1.8306762) q[1];
sx q[1];
rz(-1.13729) q[1];
rz(-1.4914289) q[3];
sx q[3];
rz(-0.36874394) q[3];
sx q[3];
rz(-0.4603995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.54625154) q[2];
sx q[2];
rz(-0.83117861) q[2];
sx q[2];
rz(1.1031995) q[2];
rz(-0.57175076) q[3];
sx q[3];
rz(-0.56643707) q[3];
sx q[3];
rz(1.5248689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6095846) q[0];
sx q[0];
rz(-0.12140618) q[0];
sx q[0];
rz(2.4914361) q[0];
rz(-1.0981759) q[1];
sx q[1];
rz(-1.2547341) q[1];
sx q[1];
rz(3.0756899) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9377334) q[0];
sx q[0];
rz(-2.2324076) q[0];
sx q[0];
rz(2.1801722) q[0];
rz(-1.6292105) q[2];
sx q[2];
rz(-1.7605392) q[2];
sx q[2];
rz(-0.74917114) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.84012953) q[1];
sx q[1];
rz(-2.3345229) q[1];
sx q[1];
rz(-2.927344) q[1];
rz(-2.8206943) q[3];
sx q[3];
rz(-0.68701744) q[3];
sx q[3];
rz(1.7248169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.99964511) q[2];
sx q[2];
rz(-1.9530752) q[2];
sx q[2];
rz(2.578242) q[2];
rz(-2.3422286) q[3];
sx q[3];
rz(-1.0146217) q[3];
sx q[3];
rz(0.25580251) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.946741) q[0];
sx q[0];
rz(-0.020337157) q[0];
sx q[0];
rz(0.49816966) q[0];
rz(0.26214552) q[1];
sx q[1];
rz(-1.405895) q[1];
sx q[1];
rz(-1.3462876) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5617843) q[0];
sx q[0];
rz(-2.0286848) q[0];
sx q[0];
rz(0.77197986) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5813229) q[2];
sx q[2];
rz(-1.5731647) q[2];
sx q[2];
rz(2.6244768) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7340103) q[1];
sx q[1];
rz(-2.2395113) q[1];
sx q[1];
rz(1.1295425) q[1];
x q[2];
rz(1.7896762) q[3];
sx q[3];
rz(-0.8026826) q[3];
sx q[3];
rz(1.7520394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4761618) q[2];
sx q[2];
rz(-2.4527145) q[2];
sx q[2];
rz(-2.2594182) q[2];
rz(-0.65794182) q[3];
sx q[3];
rz(-2.0177757) q[3];
sx q[3];
rz(-0.9995802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5049725) q[0];
sx q[0];
rz(-1.5047147) q[0];
sx q[0];
rz(2.104105) q[0];
rz(0.50719117) q[1];
sx q[1];
rz(-0.7829983) q[1];
sx q[1];
rz(-1.0601039) q[1];
rz(-1.7569594) q[2];
sx q[2];
rz(-1.2825479) q[2];
sx q[2];
rz(-2.1946236) q[2];
rz(-2.4681881) q[3];
sx q[3];
rz(-2.2670204) q[3];
sx q[3];
rz(-1.9011998) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
