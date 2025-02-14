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
rz(0.23213586) q[0];
sx q[0];
rz(-0.27096662) q[0];
sx q[0];
rz(1.969307) q[0];
rz(-1.6074033) q[1];
sx q[1];
rz(-1.489403) q[1];
sx q[1];
rz(-0.54816562) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1121824) q[0];
sx q[0];
rz(-1.8444841) q[0];
sx q[0];
rz(-0.11464439) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0127958) q[2];
sx q[2];
rz(-2.0145973) q[2];
sx q[2];
rz(1.2324126) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1480268) q[1];
sx q[1];
rz(-1.8746334) q[1];
sx q[1];
rz(0.5263473) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0788623) q[3];
sx q[3];
rz(-1.1219121) q[3];
sx q[3];
rz(-0.99454885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9008909) q[2];
sx q[2];
rz(-1.9343932) q[2];
sx q[2];
rz(1.5473676) q[2];
rz(2.2089925) q[3];
sx q[3];
rz(-2.2280732) q[3];
sx q[3];
rz(-0.80476052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43657434) q[0];
sx q[0];
rz(-2.5387634) q[0];
sx q[0];
rz(-0.21827179) q[0];
rz(-1.6328579) q[1];
sx q[1];
rz(-0.873133) q[1];
sx q[1];
rz(1.1223209) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4808761) q[0];
sx q[0];
rz(-0.97597741) q[0];
sx q[0];
rz(-0.52158458) q[0];
rz(-pi) q[1];
rz(-0.74064765) q[2];
sx q[2];
rz(-1.7193931) q[2];
sx q[2];
rz(2.1987178) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75016595) q[1];
sx q[1];
rz(-2.0199594) q[1];
sx q[1];
rz(2.5264747) q[1];
rz(-pi) q[2];
rz(1.2551941) q[3];
sx q[3];
rz(-1.3282933) q[3];
sx q[3];
rz(-0.24569233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0108769) q[2];
sx q[2];
rz(-1.9515832) q[2];
sx q[2];
rz(2.7575764) q[2];
rz(2.5859517) q[3];
sx q[3];
rz(-0.58755392) q[3];
sx q[3];
rz(-0.89239341) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4246849) q[0];
sx q[0];
rz(-0.45775828) q[0];
sx q[0];
rz(2.6440788) q[0];
rz(0.48490191) q[1];
sx q[1];
rz(-1.5496016) q[1];
sx q[1];
rz(2.7000694) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6832774) q[0];
sx q[0];
rz(-0.9181058) q[0];
sx q[0];
rz(2.9291332) q[0];
rz(-pi) q[1];
rz(0.21887987) q[2];
sx q[2];
rz(-2.0538123) q[2];
sx q[2];
rz(-2.283606) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5680629) q[1];
sx q[1];
rz(-2.2383949) q[1];
sx q[1];
rz(-0.94227643) q[1];
rz(-1.888959) q[3];
sx q[3];
rz(-1.3103879) q[3];
sx q[3];
rz(-1.3459149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.32450822) q[2];
sx q[2];
rz(-2.0980947) q[2];
sx q[2];
rz(-0.95477742) q[2];
rz(-0.020141007) q[3];
sx q[3];
rz(-0.79831278) q[3];
sx q[3];
rz(-0.32625833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37498736) q[0];
sx q[0];
rz(-2.6191481) q[0];
sx q[0];
rz(-0.0047542714) q[0];
rz(0.44113723) q[1];
sx q[1];
rz(-3.0405567) q[1];
sx q[1];
rz(-1.7316679) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72027457) q[0];
sx q[0];
rz(-1.4393115) q[0];
sx q[0];
rz(1.2340896) q[0];
rz(-0.84216046) q[2];
sx q[2];
rz(-1.7787063) q[2];
sx q[2];
rz(2.4198857) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9121352) q[1];
sx q[1];
rz(-1.3452824) q[1];
sx q[1];
rz(-2.0769801) q[1];
rz(-1.3821272) q[3];
sx q[3];
rz(-1.5975286) q[3];
sx q[3];
rz(-2.5875497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8498174) q[2];
sx q[2];
rz(-1.3128277) q[2];
sx q[2];
rz(-0.69668359) q[2];
rz(1.6126532) q[3];
sx q[3];
rz(-2.5051675) q[3];
sx q[3];
rz(2.4618885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9686389) q[0];
sx q[0];
rz(-0.30690673) q[0];
sx q[0];
rz(-0.90721834) q[0];
rz(3.0061159) q[1];
sx q[1];
rz(-1.8699346) q[1];
sx q[1];
rz(-2.4438593) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.557249) q[0];
sx q[0];
rz(-1.1797292) q[0];
sx q[0];
rz(-3.0013404) q[0];
x q[1];
rz(0.0022129755) q[2];
sx q[2];
rz(-1.7656823) q[2];
sx q[2];
rz(-2.8532078) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4042774) q[1];
sx q[1];
rz(-0.37397803) q[1];
sx q[1];
rz(-0.48980063) q[1];
rz(-pi) q[2];
rz(1.6167363) q[3];
sx q[3];
rz(-1.6005064) q[3];
sx q[3];
rz(-1.3256095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.1873143) q[2];
sx q[2];
rz(-0.66437393) q[2];
sx q[2];
rz(2.7962255) q[2];
rz(2.6464388) q[3];
sx q[3];
rz(-2.5457355) q[3];
sx q[3];
rz(-3.0143747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.0536163) q[0];
sx q[0];
rz(-1.6642267) q[0];
sx q[0];
rz(-1.9867058) q[0];
rz(-2.3467973) q[1];
sx q[1];
rz(-1.844901) q[1];
sx q[1];
rz(-2.6577139) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.799494) q[0];
sx q[0];
rz(-2.3364107) q[0];
sx q[0];
rz(2.1495887) q[0];
x q[1];
rz(2.9824798) q[2];
sx q[2];
rz(-1.0793574) q[2];
sx q[2];
rz(0.05201498) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9660873) q[1];
sx q[1];
rz(-2.6350807) q[1];
sx q[1];
rz(-1.1706074) q[1];
rz(-1.163614) q[3];
sx q[3];
rz(-2.1107499) q[3];
sx q[3];
rz(-2.0288717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.38857073) q[2];
sx q[2];
rz(-1.9066255) q[2];
sx q[2];
rz(-1.1080326) q[2];
rz(-1.5293416) q[3];
sx q[3];
rz(-3.022091) q[3];
sx q[3];
rz(2.748238) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18688467) q[0];
sx q[0];
rz(-1.0857546) q[0];
sx q[0];
rz(2.4672274) q[0];
rz(-0.87633324) q[1];
sx q[1];
rz(-2.3969711) q[1];
sx q[1];
rz(-3.1092627) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3604592) q[0];
sx q[0];
rz(-0.18258376) q[0];
sx q[0];
rz(0.83821787) q[0];
rz(-0.22458073) q[2];
sx q[2];
rz(-1.7214643) q[2];
sx q[2];
rz(1.2008787) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9385953) q[1];
sx q[1];
rz(-1.055966) q[1];
sx q[1];
rz(1.1802303) q[1];
rz(-pi) q[2];
rz(2.1635734) q[3];
sx q[3];
rz(-1.4109352) q[3];
sx q[3];
rz(0.39409742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9548367) q[2];
sx q[2];
rz(-2.7080471) q[2];
sx q[2];
rz(-1.0369066) q[2];
rz(1.9131276) q[3];
sx q[3];
rz(-0.90932536) q[3];
sx q[3];
rz(-0.19074805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4036338) q[0];
sx q[0];
rz(-3.1138595) q[0];
sx q[0];
rz(0.24895689) q[0];
rz(0.20949334) q[1];
sx q[1];
rz(-1.341235) q[1];
sx q[1];
rz(0.6627717) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4512206) q[0];
sx q[0];
rz(-1.5826384) q[0];
sx q[0];
rz(-1.5607587) q[0];
x q[1];
rz(-2.6423062) q[2];
sx q[2];
rz(-2.3134276) q[2];
sx q[2];
rz(1.7940831) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.31964916) q[1];
sx q[1];
rz(-1.3256097) q[1];
sx q[1];
rz(2.3169434) q[1];
rz(-0.39318496) q[3];
sx q[3];
rz(-0.69614702) q[3];
sx q[3];
rz(-0.37567155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2304307) q[2];
sx q[2];
rz(-1.570236) q[2];
sx q[2];
rz(-2.7479808) q[2];
rz(-0.39129928) q[3];
sx q[3];
rz(-2.6063882) q[3];
sx q[3];
rz(2.3894501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0126208) q[0];
sx q[0];
rz(-0.18652815) q[0];
sx q[0];
rz(2.5616126) q[0];
rz(-2.773556) q[1];
sx q[1];
rz(-2.3263704) q[1];
sx q[1];
rz(-0.11485242) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2987679) q[0];
sx q[0];
rz(-1.4359763) q[0];
sx q[0];
rz(2.299286) q[0];
x q[1];
rz(1.2602368) q[2];
sx q[2];
rz(-2.680696) q[2];
sx q[2];
rz(0.68379921) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8424994) q[1];
sx q[1];
rz(-2.7242492) q[1];
sx q[1];
rz(2.4853766) q[1];
rz(0.87646342) q[3];
sx q[3];
rz(-1.8871347) q[3];
sx q[3];
rz(-2.9522459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4174691) q[2];
sx q[2];
rz(-2.8026411) q[2];
sx q[2];
rz(0.69619703) q[2];
rz(-1.1160858) q[3];
sx q[3];
rz(-2.3510272) q[3];
sx q[3];
rz(-0.74151403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8490132) q[0];
sx q[0];
rz(-0.96939033) q[0];
sx q[0];
rz(-0.27266362) q[0];
rz(-2.4478681) q[1];
sx q[1];
rz(-1.1169746) q[1];
sx q[1];
rz(-0.56347096) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7359223) q[0];
sx q[0];
rz(-1.6077245) q[0];
sx q[0];
rz(-1.4775299) q[0];
rz(-pi) q[1];
rz(2.5921037) q[2];
sx q[2];
rz(-2.8656883) q[2];
sx q[2];
rz(-0.23978182) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2595926) q[1];
sx q[1];
rz(-1.2597076) q[1];
sx q[1];
rz(3.056219) q[1];
rz(-pi) q[2];
rz(0.93592398) q[3];
sx q[3];
rz(-2.1800123) q[3];
sx q[3];
rz(0.09906957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1571265) q[2];
sx q[2];
rz(-1.0447341) q[2];
sx q[2];
rz(-0.54894471) q[2];
rz(-0.99556154) q[3];
sx q[3];
rz(-0.82058161) q[3];
sx q[3];
rz(0.63001776) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7231049) q[0];
sx q[0];
rz(-0.99011078) q[0];
sx q[0];
rz(-0.44575442) q[0];
rz(2.2743629) q[1];
sx q[1];
rz(-1.1031311) q[1];
sx q[1];
rz(1.1581609) q[1];
rz(0.42648496) q[2];
sx q[2];
rz(-2.6250962) q[2];
sx q[2];
rz(0.442183) q[2];
rz(1.7165202) q[3];
sx q[3];
rz(-1.0276262) q[3];
sx q[3];
rz(0.75158391) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
