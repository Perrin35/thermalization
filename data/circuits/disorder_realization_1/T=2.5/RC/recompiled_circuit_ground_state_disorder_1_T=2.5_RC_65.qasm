OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.98439944) q[0];
sx q[0];
rz(-1.1097044) q[0];
sx q[0];
rz(-2.1362526) q[0];
rz(-2.4069064) q[1];
sx q[1];
rz(-0.35597304) q[1];
sx q[1];
rz(-2.2495143) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1482502) q[0];
sx q[0];
rz(-2.3138232) q[0];
sx q[0];
rz(-1.7048852) q[0];
rz(-pi) q[1];
rz(2.1744035) q[2];
sx q[2];
rz(-1.693486) q[2];
sx q[2];
rz(1.762378) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5106543) q[1];
sx q[1];
rz(-1.0882241) q[1];
sx q[1];
rz(1.4605182) q[1];
rz(-pi) q[2];
rz(-2.0964811) q[3];
sx q[3];
rz(-2.1108004) q[3];
sx q[3];
rz(-1.2837376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4165108) q[2];
sx q[2];
rz(-1.9981013) q[2];
sx q[2];
rz(1.4408646) q[2];
rz(-0.95669389) q[3];
sx q[3];
rz(-1.1172833) q[3];
sx q[3];
rz(0.24334894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44593921) q[0];
sx q[0];
rz(-2.74701) q[0];
sx q[0];
rz(1.12895) q[0];
rz(0.24540643) q[1];
sx q[1];
rz(-1.8281728) q[1];
sx q[1];
rz(-2.479877) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.322987) q[0];
sx q[0];
rz(-1.855369) q[0];
sx q[0];
rz(-2.5045082) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3459588) q[2];
sx q[2];
rz(-1.1974466) q[2];
sx q[2];
rz(2.7527347) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0894818) q[1];
sx q[1];
rz(-2.8683337) q[1];
sx q[1];
rz(0.11788003) q[1];
rz(-pi) q[2];
rz(-1.861694) q[3];
sx q[3];
rz(-1.6637633) q[3];
sx q[3];
rz(2.7806843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5976065) q[2];
sx q[2];
rz(-0.49763766) q[2];
sx q[2];
rz(1.0428766) q[2];
rz(1.6889702) q[3];
sx q[3];
rz(-0.85113227) q[3];
sx q[3];
rz(2.0378621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3747028) q[0];
sx q[0];
rz(-0.68793982) q[0];
sx q[0];
rz(-1.3091298) q[0];
rz(0.4492999) q[1];
sx q[1];
rz(-2.4520912) q[1];
sx q[1];
rz(1.4264533) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2319205) q[0];
sx q[0];
rz(-1.3065152) q[0];
sx q[0];
rz(2.8876165) q[0];
rz(-3.1235891) q[2];
sx q[2];
rz(-1.1483542) q[2];
sx q[2];
rz(-1.9033592) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.92628463) q[1];
sx q[1];
rz(-2.4325772) q[1];
sx q[1];
rz(-0.11911094) q[1];
x q[2];
rz(-3.0578311) q[3];
sx q[3];
rz(-1.8159831) q[3];
sx q[3];
rz(0.29473588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7766777) q[2];
sx q[2];
rz(-1.1740351) q[2];
sx q[2];
rz(-3.0752227) q[2];
rz(2.5868609) q[3];
sx q[3];
rz(-1.4812255) q[3];
sx q[3];
rz(1.6248645) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9471112) q[0];
sx q[0];
rz(-0.28488657) q[0];
sx q[0];
rz(2.3696005) q[0];
rz(-0.66328612) q[1];
sx q[1];
rz(-0.74378496) q[1];
sx q[1];
rz(-3.0139121) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4366764) q[0];
sx q[0];
rz(-0.66670115) q[0];
sx q[0];
rz(0.0078860869) q[0];
rz(-pi) q[1];
rz(-0.56945412) q[2];
sx q[2];
rz(-0.63210154) q[2];
sx q[2];
rz(-2.5394627) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7230715) q[1];
sx q[1];
rz(-1.1643895) q[1];
sx q[1];
rz(-1.0999098) q[1];
rz(-1.4127619) q[3];
sx q[3];
rz(-1.1452066) q[3];
sx q[3];
rz(1.1273718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6105303) q[2];
sx q[2];
rz(-2.1805111) q[2];
sx q[2];
rz(-2.5430211) q[2];
rz(-0.5851723) q[3];
sx q[3];
rz(-1.7681311) q[3];
sx q[3];
rz(-3.0857871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7912306) q[0];
sx q[0];
rz(-1.5086011) q[0];
sx q[0];
rz(-0.97306657) q[0];
rz(0.19634253) q[1];
sx q[1];
rz(-2.0025496) q[1];
sx q[1];
rz(-1.4686718) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3298686) q[0];
sx q[0];
rz(-1.1783333) q[0];
sx q[0];
rz(-2.3377665) q[0];
x q[1];
rz(1.0650915) q[2];
sx q[2];
rz(-1.347216) q[2];
sx q[2];
rz(-2.2044942) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0397075) q[1];
sx q[1];
rz(-0.3488003) q[1];
sx q[1];
rz(1.0451911) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0319388) q[3];
sx q[3];
rz(-1.6183637) q[3];
sx q[3];
rz(2.4470234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.68957442) q[2];
sx q[2];
rz(-1.8147899) q[2];
sx q[2];
rz(-0.16239521) q[2];
rz(1.1299805) q[3];
sx q[3];
rz(-0.88310784) q[3];
sx q[3];
rz(-1.1348772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3910386) q[0];
sx q[0];
rz(-2.6226608) q[0];
sx q[0];
rz(-3.0911875) q[0];
rz(-2.9734036) q[1];
sx q[1];
rz(-1.5610118) q[1];
sx q[1];
rz(0.87356299) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03913232) q[0];
sx q[0];
rz(-0.96193571) q[0];
sx q[0];
rz(1.8692506) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1077439) q[2];
sx q[2];
rz(-1.7122972) q[2];
sx q[2];
rz(-1.3901097) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7616523) q[1];
sx q[1];
rz(-0.72387513) q[1];
sx q[1];
rz(0.70751247) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5477375) q[3];
sx q[3];
rz(-2.7866413) q[3];
sx q[3];
rz(1.5221727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0198589) q[2];
sx q[2];
rz(-2.9254318) q[2];
sx q[2];
rz(-2.7609694) q[2];
rz(-2.5753283) q[3];
sx q[3];
rz(-1.7411313) q[3];
sx q[3];
rz(2.3316021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5009163) q[0];
sx q[0];
rz(-0.2114507) q[0];
sx q[0];
rz(2.7389615) q[0];
rz(-1.3817878) q[1];
sx q[1];
rz(-2.333162) q[1];
sx q[1];
rz(-3.0852539) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9531247) q[0];
sx q[0];
rz(-1.2446277) q[0];
sx q[0];
rz(-1.1475546) q[0];
rz(-pi) q[1];
rz(-0.041002657) q[2];
sx q[2];
rz(-2.0411388) q[2];
sx q[2];
rz(-0.13861632) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.90625396) q[1];
sx q[1];
rz(-1.5829931) q[1];
sx q[1];
rz(-1.0066628) q[1];
rz(-2.4640637) q[3];
sx q[3];
rz(-1.9871681) q[3];
sx q[3];
rz(0.29034055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.23400433) q[2];
sx q[2];
rz(-1.3434709) q[2];
sx q[2];
rz(2.6410356) q[2];
rz(3.0808926) q[3];
sx q[3];
rz(-1.3541636) q[3];
sx q[3];
rz(2.6589656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.56977) q[0];
sx q[0];
rz(-1.7800542) q[0];
sx q[0];
rz(1.7806336) q[0];
rz(0.743615) q[1];
sx q[1];
rz(-1.7230325) q[1];
sx q[1];
rz(3.0027711) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1847938) q[0];
sx q[0];
rz(-1.0652055) q[0];
sx q[0];
rz(0.57539815) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1651504) q[2];
sx q[2];
rz(-1.1757554) q[2];
sx q[2];
rz(0.76752418) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.239153) q[1];
sx q[1];
rz(-0.3867074) q[1];
sx q[1];
rz(1.0104806) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0431792) q[3];
sx q[3];
rz(-2.7076027) q[3];
sx q[3];
rz(1.4066309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4678141) q[2];
sx q[2];
rz(-2.0970586) q[2];
sx q[2];
rz(-0.71411258) q[2];
rz(2.1821187) q[3];
sx q[3];
rz(-1.6474612) q[3];
sx q[3];
rz(-0.97904557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4176843) q[0];
sx q[0];
rz(-2.4651616) q[0];
sx q[0];
rz(-3.0133001) q[0];
rz(-2.7543606) q[1];
sx q[1];
rz(-2.5435244) q[1];
sx q[1];
rz(-1.3234352) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0664862) q[0];
sx q[0];
rz(-1.3818701) q[0];
sx q[0];
rz(-0.099204258) q[0];
rz(-pi) q[1];
rz(2.2567184) q[2];
sx q[2];
rz(-1.0909683) q[2];
sx q[2];
rz(-0.72140933) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.379017) q[1];
sx q[1];
rz(-1.5847995) q[1];
sx q[1];
rz(1.4268616) q[1];
x q[2];
rz(-0.22289688) q[3];
sx q[3];
rz(-0.74046293) q[3];
sx q[3];
rz(-2.2929913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6748176) q[2];
sx q[2];
rz(-1.520949) q[2];
sx q[2];
rz(-2.9368994) q[2];
rz(-1.8681059) q[3];
sx q[3];
rz(-0.73015648) q[3];
sx q[3];
rz(1.5761121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8650763) q[0];
sx q[0];
rz(-0.62224329) q[0];
sx q[0];
rz(-2.9283071) q[0];
rz(-2.9215096) q[1];
sx q[1];
rz(-1.8890231) q[1];
sx q[1];
rz(1.3594886) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8187788) q[0];
sx q[0];
rz(-1.5867763) q[0];
sx q[0];
rz(0.009601618) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47410902) q[2];
sx q[2];
rz(-1.7987937) q[2];
sx q[2];
rz(2.1257509) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9656187) q[1];
sx q[1];
rz(-0.63116628) q[1];
sx q[1];
rz(0.96221535) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1207934) q[3];
sx q[3];
rz(-0.74899835) q[3];
sx q[3];
rz(1.7089546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6243817) q[2];
sx q[2];
rz(-1.5959847) q[2];
sx q[2];
rz(0.33779302) q[2];
rz(1.7289915) q[3];
sx q[3];
rz(-1.1510886) q[3];
sx q[3];
rz(-0.82144773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74078858) q[0];
sx q[0];
rz(-2.3176226) q[0];
sx q[0];
rz(1.2299706) q[0];
rz(-0.38372718) q[1];
sx q[1];
rz(-1.6576672) q[1];
sx q[1];
rz(-2.3794649) q[1];
rz(2.7376851) q[2];
sx q[2];
rz(-1.8941034) q[2];
sx q[2];
rz(-0.20878172) q[2];
rz(-1.0596342) q[3];
sx q[3];
rz(-1.3564902) q[3];
sx q[3];
rz(-1.0839957) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
