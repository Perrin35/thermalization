OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8813397) q[0];
sx q[0];
rz(-0.94085675) q[0];
sx q[0];
rz(2.9139304) q[0];
rz(0.28336278) q[1];
sx q[1];
rz(3.5609666) q[1];
sx q[1];
rz(11.279439) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44449466) q[0];
sx q[0];
rz(-2.6711406) q[0];
sx q[0];
rz(1.6906428) q[0];
rz(-pi) q[1];
rz(0.46704328) q[2];
sx q[2];
rz(-1.9775632) q[2];
sx q[2];
rz(1.1554935) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70482774) q[1];
sx q[1];
rz(-0.65751644) q[1];
sx q[1];
rz(-2.2566811) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8836796) q[3];
sx q[3];
rz(-0.95525026) q[3];
sx q[3];
rz(-0.93759495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2490354) q[2];
sx q[2];
rz(-2.1302569) q[2];
sx q[2];
rz(2.2300143) q[2];
rz(2.3859731) q[3];
sx q[3];
rz(-2.8204462) q[3];
sx q[3];
rz(2.6452981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(2.3590473) q[0];
sx q[0];
rz(-1.8308715) q[0];
sx q[0];
rz(2.0027335) q[0];
rz(0.62659872) q[1];
sx q[1];
rz(-2.7732924) q[1];
sx q[1];
rz(2.0638594) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2127578) q[0];
sx q[0];
rz(-0.99883119) q[0];
sx q[0];
rz(-1.1642745) q[0];
x q[1];
rz(-1.1280649) q[2];
sx q[2];
rz(-1.2473543) q[2];
sx q[2];
rz(-1.6328822) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6328395) q[1];
sx q[1];
rz(-1.6378073) q[1];
sx q[1];
rz(-2.8772023) q[1];
x q[2];
rz(2.0502362) q[3];
sx q[3];
rz(-0.78244996) q[3];
sx q[3];
rz(-1.2408011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.44800147) q[2];
sx q[2];
rz(-2.3841264) q[2];
sx q[2];
rz(-0.21270154) q[2];
rz(1.5564144) q[3];
sx q[3];
rz(-0.74376619) q[3];
sx q[3];
rz(1.0630382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14146516) q[0];
sx q[0];
rz(-0.26697049) q[0];
sx q[0];
rz(2.2947327) q[0];
rz(-0.2521387) q[1];
sx q[1];
rz(-0.82770258) q[1];
sx q[1];
rz(-2.491378) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0747809) q[0];
sx q[0];
rz(-1.909314) q[0];
sx q[0];
rz(1.6174497) q[0];
rz(0.89805897) q[2];
sx q[2];
rz(-0.89623797) q[2];
sx q[2];
rz(-2.7314699) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0155772) q[1];
sx q[1];
rz(-1.6815876) q[1];
sx q[1];
rz(1.2726102) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7103151) q[3];
sx q[3];
rz(-0.85429885) q[3];
sx q[3];
rz(2.7950479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3113159) q[2];
sx q[2];
rz(-2.4429784) q[2];
sx q[2];
rz(-0.96381956) q[2];
rz(-1.7802995) q[3];
sx q[3];
rz(-0.45069525) q[3];
sx q[3];
rz(3.0986339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90218246) q[0];
sx q[0];
rz(-2.9189411) q[0];
sx q[0];
rz(-2.1233001) q[0];
rz(-2.2564383) q[1];
sx q[1];
rz(-2.8926909) q[1];
sx q[1];
rz(-2.8525066) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1219629) q[0];
sx q[0];
rz(-1.2696854) q[0];
sx q[0];
rz(2.798978) q[0];
rz(-1.1491597) q[2];
sx q[2];
rz(-0.9767864) q[2];
sx q[2];
rz(-0.27911148) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8377467) q[1];
sx q[1];
rz(-1.0582542) q[1];
sx q[1];
rz(-0.20767494) q[1];
x q[2];
rz(0.21896514) q[3];
sx q[3];
rz(-1.4462399) q[3];
sx q[3];
rz(-1.3606461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.79124147) q[2];
sx q[2];
rz(-1.505932) q[2];
sx q[2];
rz(-1.2887456) q[2];
rz(-0.57981235) q[3];
sx q[3];
rz(-2.6668187) q[3];
sx q[3];
rz(0.60540664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8346005) q[0];
sx q[0];
rz(-0.51486105) q[0];
sx q[0];
rz(-2.8173764) q[0];
rz(3.1027555) q[1];
sx q[1];
rz(-2.4034998) q[1];
sx q[1];
rz(-2.0447581) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0457458) q[0];
sx q[0];
rz(-2.7660094) q[0];
sx q[0];
rz(-2.4235241) q[0];
rz(-pi) q[1];
rz(1.4154424) q[2];
sx q[2];
rz(-0.58621472) q[2];
sx q[2];
rz(-1.6833351) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8406163) q[1];
sx q[1];
rz(-2.8368717) q[1];
sx q[1];
rz(-1.4304377) q[1];
rz(-1.3665133) q[3];
sx q[3];
rz(-0.42821233) q[3];
sx q[3];
rz(-2.7985559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0754806) q[2];
sx q[2];
rz(-0.26618633) q[2];
sx q[2];
rz(-3.10293) q[2];
rz(-1.4085116) q[3];
sx q[3];
rz(-1.446529) q[3];
sx q[3];
rz(-0.2814289) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51327813) q[0];
sx q[0];
rz(-0.81858855) q[0];
sx q[0];
rz(1.0191089) q[0];
rz(-0.45711532) q[1];
sx q[1];
rz(-2.8725084) q[1];
sx q[1];
rz(-1.7105182) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015648677) q[0];
sx q[0];
rz(-1.5971184) q[0];
sx q[0];
rz(0.20080815) q[0];
rz(1.9408631) q[2];
sx q[2];
rz(-2.5197755) q[2];
sx q[2];
rz(-1.3696826) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6999457) q[1];
sx q[1];
rz(-1.4829548) q[1];
sx q[1];
rz(-1.2097174) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7442211) q[3];
sx q[3];
rz(-0.89029593) q[3];
sx q[3];
rz(1.3272663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.21923253) q[2];
sx q[2];
rz(-1.7722426) q[2];
sx q[2];
rz(-2.5617981) q[2];
rz(-0.29948768) q[3];
sx q[3];
rz(-2.6087285) q[3];
sx q[3];
rz(-1.2286435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86094588) q[0];
sx q[0];
rz(-1.6609284) q[0];
sx q[0];
rz(3.0378367) q[0];
rz(-2.5170028) q[1];
sx q[1];
rz(-0.92085212) q[1];
sx q[1];
rz(-1.7594899) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9860024) q[0];
sx q[0];
rz(-2.3891923) q[0];
sx q[0];
rz(-1.2563906) q[0];
rz(-pi) q[1];
rz(1.4598249) q[2];
sx q[2];
rz(-1.2726415) q[2];
sx q[2];
rz(-2.0046749) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4098072) q[1];
sx q[1];
rz(-2.5587132) q[1];
sx q[1];
rz(-1.0686841) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2861038) q[3];
sx q[3];
rz(-2.2119388) q[3];
sx q[3];
rz(-0.59635983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.78911191) q[2];
sx q[2];
rz(-1.666297) q[2];
sx q[2];
rz(0.86461198) q[2];
rz(2.9027446) q[3];
sx q[3];
rz(-0.75967234) q[3];
sx q[3];
rz(-2.4281832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9637941) q[0];
sx q[0];
rz(-0.56121427) q[0];
sx q[0];
rz(-2.906565) q[0];
rz(0.66559732) q[1];
sx q[1];
rz(-1.1382297) q[1];
sx q[1];
rz(-1.4733018) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4427933) q[0];
sx q[0];
rz(-0.59663749) q[0];
sx q[0];
rz(1.9477316) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8741415) q[2];
sx q[2];
rz(-0.48201928) q[2];
sx q[2];
rz(-2.6798525) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7104147) q[1];
sx q[1];
rz(-1.3873552) q[1];
sx q[1];
rz(1.6448433) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5258011) q[3];
sx q[3];
rz(-2.0215086) q[3];
sx q[3];
rz(-1.352528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.16234806) q[2];
sx q[2];
rz(-2.9408216) q[2];
sx q[2];
rz(-0.46094224) q[2];
rz(1.0848684) q[3];
sx q[3];
rz(-2.3960787) q[3];
sx q[3];
rz(-2.357024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(0.064086176) q[0];
sx q[0];
rz(-2.6106847) q[0];
sx q[0];
rz(-2.532646) q[0];
rz(-0.56600904) q[1];
sx q[1];
rz(-1.1606263) q[1];
sx q[1];
rz(1.2014679) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5574359) q[0];
sx q[0];
rz(-0.13579255) q[0];
sx q[0];
rz(-1.5103755) q[0];
x q[1];
rz(2.9329691) q[2];
sx q[2];
rz(-2.0498599) q[2];
sx q[2];
rz(2.5076659) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4200099) q[1];
sx q[1];
rz(-2.3370565) q[1];
sx q[1];
rz(-1.700042) q[1];
rz(0.5769095) q[3];
sx q[3];
rz(-0.52191496) q[3];
sx q[3];
rz(-1.0229223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.10298771) q[2];
sx q[2];
rz(-1.0706341) q[2];
sx q[2];
rz(-1.1919682) q[2];
rz(2.1206756) q[3];
sx q[3];
rz(-0.20191419) q[3];
sx q[3];
rz(-1.5405704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22333071) q[0];
sx q[0];
rz(-2.9409565) q[0];
sx q[0];
rz(-0.9675135) q[0];
rz(-1.4406904) q[1];
sx q[1];
rz(-0.43165019) q[1];
sx q[1];
rz(2.7593625) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0959045) q[0];
sx q[0];
rz(-1.75215) q[0];
sx q[0];
rz(-2.1482094) q[0];
rz(-pi) q[1];
rz(-0.43469825) q[2];
sx q[2];
rz(-1.189173) q[2];
sx q[2];
rz(-1.3023072) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6401419) q[1];
sx q[1];
rz(-0.79265187) q[1];
sx q[1];
rz(-0.16736302) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54602539) q[3];
sx q[3];
rz(-1.556753) q[3];
sx q[3];
rz(1.394681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.091577) q[2];
sx q[2];
rz(-2.2735368) q[2];
sx q[2];
rz(-0.78563219) q[2];
rz(0.026963726) q[3];
sx q[3];
rz(-1.5188768) q[3];
sx q[3];
rz(-0.54076076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3383269) q[0];
sx q[0];
rz(-1.6852408) q[0];
sx q[0];
rz(2.2483873) q[0];
rz(0.66435736) q[1];
sx q[1];
rz(-1.9021481) q[1];
sx q[1];
rz(-0.93217168) q[1];
rz(0.33334941) q[2];
sx q[2];
rz(-0.7569505) q[2];
sx q[2];
rz(-1.0406582) q[2];
rz(-2.5369) q[3];
sx q[3];
rz(-1.562955) q[3];
sx q[3];
rz(2.4969586) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
