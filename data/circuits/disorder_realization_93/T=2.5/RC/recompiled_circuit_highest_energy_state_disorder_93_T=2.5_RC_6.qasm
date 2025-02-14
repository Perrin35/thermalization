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
rz(0.33557284) q[0];
sx q[0];
rz(-1.4541452) q[0];
sx q[0];
rz(-2.1079221) q[0];
rz(4.5728788) q[1];
sx q[1];
rz(3.6990777) q[1];
sx q[1];
rz(7.2929444) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5752707) q[0];
sx q[0];
rz(-0.22804865) q[0];
sx q[0];
rz(2.1581677) q[0];
rz(1.4883243) q[2];
sx q[2];
rz(-1.5670652) q[2];
sx q[2];
rz(-0.53860215) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3897455) q[1];
sx q[1];
rz(-1.6865786) q[1];
sx q[1];
rz(0.16645653) q[1];
x q[2];
rz(-1.5663773) q[3];
sx q[3];
rz(-1.7130245) q[3];
sx q[3];
rz(-2.9517236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7734163) q[2];
sx q[2];
rz(-2.2645617) q[2];
sx q[2];
rz(0.82484335) q[2];
rz(-2.8090737) q[3];
sx q[3];
rz(-2.2907292) q[3];
sx q[3];
rz(2.6578145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2291718) q[0];
sx q[0];
rz(-0.084675463) q[0];
sx q[0];
rz(0.53572768) q[0];
rz(2.3748705) q[1];
sx q[1];
rz(-1.7886536) q[1];
sx q[1];
rz(-1.3923233) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.509393) q[0];
sx q[0];
rz(-2.0344816) q[0];
sx q[0];
rz(0.38393387) q[0];
rz(1.5803945) q[2];
sx q[2];
rz(-1.9119153) q[2];
sx q[2];
rz(0.3915325) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.630188) q[1];
sx q[1];
rz(-1.5923481) q[1];
sx q[1];
rz(0.81540108) q[1];
rz(-1.6757319) q[3];
sx q[3];
rz(-1.6989048) q[3];
sx q[3];
rz(-0.48612938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.9048189) q[2];
sx q[2];
rz(-0.83319131) q[2];
sx q[2];
rz(1.0986249) q[2];
rz(2.3084579) q[3];
sx q[3];
rz(-1.5980709) q[3];
sx q[3];
rz(-1.0540849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58040923) q[0];
sx q[0];
rz(-1.1495178) q[0];
sx q[0];
rz(-2.4533601) q[0];
rz(-1.2511823) q[1];
sx q[1];
rz(-2.4871608) q[1];
sx q[1];
rz(1.2122663) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1935342) q[0];
sx q[0];
rz(-2.2431264) q[0];
sx q[0];
rz(-2.1904519) q[0];
rz(-pi) q[1];
rz(2.9040292) q[2];
sx q[2];
rz(-1.8389987) q[2];
sx q[2];
rz(-0.99246565) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7408091) q[1];
sx q[1];
rz(-2.2297528) q[1];
sx q[1];
rz(-1.5636921) q[1];
rz(-1.5098677) q[3];
sx q[3];
rz(-0.079515545) q[3];
sx q[3];
rz(-2.7729976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7324098) q[2];
sx q[2];
rz(-0.99486351) q[2];
sx q[2];
rz(2.691972) q[2];
rz(2.2879587) q[3];
sx q[3];
rz(-0.64695224) q[3];
sx q[3];
rz(2.4322815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7483826) q[0];
sx q[0];
rz(-1.0164096) q[0];
sx q[0];
rz(-2.7384695) q[0];
rz(2.368811) q[1];
sx q[1];
rz(-1.9085596) q[1];
sx q[1];
rz(-0.82537878) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3365509) q[0];
sx q[0];
rz(-2.2580349) q[0];
sx q[0];
rz(-2.7206793) q[0];
rz(-2.4140777) q[2];
sx q[2];
rz(-3.10559) q[2];
sx q[2];
rz(1.1628448) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9627208) q[1];
sx q[1];
rz(-2.555518) q[1];
sx q[1];
rz(-1.7078215) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.585936) q[3];
sx q[3];
rz(-2.1534277) q[3];
sx q[3];
rz(0.61911303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0394502) q[2];
sx q[2];
rz(-1.2692229) q[2];
sx q[2];
rz(-2.6329182) q[2];
rz(1.8968808) q[3];
sx q[3];
rz(-1.0735984) q[3];
sx q[3];
rz(-1.8083474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4077069) q[0];
sx q[0];
rz(-1.5835967) q[0];
sx q[0];
rz(2.7852614) q[0];
rz(-1.7286667) q[1];
sx q[1];
rz(-1.3099542) q[1];
sx q[1];
rz(1.8467356) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39065162) q[0];
sx q[0];
rz(-1.3518419) q[0];
sx q[0];
rz(0.44415565) q[0];
rz(-pi) q[1];
rz(0.21759502) q[2];
sx q[2];
rz(-2.4684973) q[2];
sx q[2];
rz(2.4270647) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9425671) q[1];
sx q[1];
rz(-0.46442214) q[1];
sx q[1];
rz(-2.6075122) q[1];
rz(-pi) q[2];
rz(0.89957063) q[3];
sx q[3];
rz(-2.482467) q[3];
sx q[3];
rz(2.5289583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9406708) q[2];
sx q[2];
rz(-1.7565497) q[2];
sx q[2];
rz(0.6400288) q[2];
rz(-0.43404964) q[3];
sx q[3];
rz(-0.95165747) q[3];
sx q[3];
rz(-2.0210338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5459179) q[0];
sx q[0];
rz(-0.43788236) q[0];
sx q[0];
rz(-0.75089279) q[0];
rz(-0.76260507) q[1];
sx q[1];
rz(-1.4354939) q[1];
sx q[1];
rz(0.5336175) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5845244) q[0];
sx q[0];
rz(-1.2831339) q[0];
sx q[0];
rz(2.8329222) q[0];
rz(-pi) q[1];
rz(2.7087174) q[2];
sx q[2];
rz(-0.85202571) q[2];
sx q[2];
rz(-1.5804039) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.74063674) q[1];
sx q[1];
rz(-1.6888685) q[1];
sx q[1];
rz(2.2406039) q[1];
rz(-pi) q[2];
rz(-2.3271766) q[3];
sx q[3];
rz(-0.73668639) q[3];
sx q[3];
rz(1.8131922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.077347191) q[2];
sx q[2];
rz(-2.7392445) q[2];
sx q[2];
rz(2.2557491) q[2];
rz(-1.3028076) q[3];
sx q[3];
rz(-1.1833444) q[3];
sx q[3];
rz(-0.47811374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.2122022) q[0];
sx q[0];
rz(-0.92340702) q[0];
sx q[0];
rz(-2.5481664) q[0];
rz(1.5133096) q[1];
sx q[1];
rz(-1.9624886) q[1];
sx q[1];
rz(2.5679307) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5712515) q[0];
sx q[0];
rz(-1.9431337) q[0];
sx q[0];
rz(-0.516173) q[0];
x q[1];
rz(-0.57777496) q[2];
sx q[2];
rz(-1.7794357) q[2];
sx q[2];
rz(1.278217) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1036914) q[1];
sx q[1];
rz(-1.173844) q[1];
sx q[1];
rz(1.3034225) q[1];
x q[2];
rz(0.078646831) q[3];
sx q[3];
rz(-1.5203195) q[3];
sx q[3];
rz(0.61928643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.71184984) q[2];
sx q[2];
rz(-1.9387127) q[2];
sx q[2];
rz(1.4329866) q[2];
rz(1.9728707) q[3];
sx q[3];
rz(-2.70372) q[3];
sx q[3];
rz(1.6247113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2699921) q[0];
sx q[0];
rz(-0.90684909) q[0];
sx q[0];
rz(-0.20371833) q[0];
rz(1.3153971) q[1];
sx q[1];
rz(-1.8351277) q[1];
sx q[1];
rz(-1.3320097) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7527134) q[0];
sx q[0];
rz(-0.95001555) q[0];
sx q[0];
rz(-2.3546773) q[0];
rz(-0.88896759) q[2];
sx q[2];
rz(-0.79304129) q[2];
sx q[2];
rz(2.5630132) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.89573025) q[1];
sx q[1];
rz(-1.9612211) q[1];
sx q[1];
rz(-0.67710442) q[1];
x q[2];
rz(0.90059728) q[3];
sx q[3];
rz(-0.20684563) q[3];
sx q[3];
rz(1.2634009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3709162) q[2];
sx q[2];
rz(-0.081826536) q[2];
sx q[2];
rz(-1.7601684) q[2];
rz(0.55111876) q[3];
sx q[3];
rz(-1.328732) q[3];
sx q[3];
rz(-0.74530017) q[3];
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
rz(1.7119174) q[0];
sx q[0];
rz(-1.6918007) q[0];
sx q[0];
rz(1.8823189) q[0];
rz(1.7970386) q[1];
sx q[1];
rz(-0.63251907) q[1];
sx q[1];
rz(-1.9154027) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38905242) q[0];
sx q[0];
rz(-1.8551833) q[0];
sx q[0];
rz(3.000598) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38649747) q[2];
sx q[2];
rz(-1.1697239) q[2];
sx q[2];
rz(2.3569466) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.74599712) q[1];
sx q[1];
rz(-1.5052553) q[1];
sx q[1];
rz(-2.3941674) q[1];
rz(2.8060447) q[3];
sx q[3];
rz(-2.0651544) q[3];
sx q[3];
rz(1.7691607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.58174497) q[2];
sx q[2];
rz(-0.86220828) q[2];
sx q[2];
rz(-1.6017412) q[2];
rz(1.7848484) q[3];
sx q[3];
rz(-0.53250766) q[3];
sx q[3];
rz(1.2342359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(2.8643879) q[0];
sx q[0];
rz(-1.5437523) q[0];
sx q[0];
rz(0.10449115) q[0];
rz(-1.3970207) q[1];
sx q[1];
rz(-1.7897768) q[1];
sx q[1];
rz(1.5207965) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8408729) q[0];
sx q[0];
rz(-1.3618293) q[0];
sx q[0];
rz(0.74700256) q[0];
rz(0.42664032) q[2];
sx q[2];
rz(-0.37084118) q[2];
sx q[2];
rz(-1.1877354) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4222153) q[1];
sx q[1];
rz(-1.4669734) q[1];
sx q[1];
rz(-2.8524998) q[1];
rz(-pi) q[2];
rz(1.4827221) q[3];
sx q[3];
rz(-2.5582426) q[3];
sx q[3];
rz(-0.47881918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.70667679) q[2];
sx q[2];
rz(-0.84182635) q[2];
sx q[2];
rz(-1.288877) q[2];
rz(-2.1364818) q[3];
sx q[3];
rz(-1.0594599) q[3];
sx q[3];
rz(3.03481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91916753) q[0];
sx q[0];
rz(-1.5821624) q[0];
sx q[0];
rz(-0.17056175) q[0];
rz(-0.87986058) q[1];
sx q[1];
rz(-0.51645551) q[1];
sx q[1];
rz(-2.5194306) q[1];
rz(0.11779412) q[2];
sx q[2];
rz(-0.5024903) q[2];
sx q[2];
rz(2.6033664) q[2];
rz(-2.5414657) q[3];
sx q[3];
rz(-1.5059581) q[3];
sx q[3];
rz(-3.0653421) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
