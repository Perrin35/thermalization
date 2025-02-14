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
rz(-1.0626592) q[0];
sx q[0];
rz(-0.81544977) q[0];
sx q[0];
rz(0.71001473) q[0];
rz(-0.31399909) q[1];
sx q[1];
rz(-0.93500885) q[1];
sx q[1];
rz(11.234565) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79373097) q[0];
sx q[0];
rz(-1.2958741) q[0];
sx q[0];
rz(1.8040276) q[0];
rz(2.7983886) q[2];
sx q[2];
rz(-1.8095952) q[2];
sx q[2];
rz(-2.3607139) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9985245) q[1];
sx q[1];
rz(-1.3245163) q[1];
sx q[1];
rz(0.95264901) q[1];
rz(2.2043742) q[3];
sx q[3];
rz(-1.9301747) q[3];
sx q[3];
rz(1.0868974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.72267246) q[2];
sx q[2];
rz(-1.924943) q[2];
sx q[2];
rz(-0.81532064) q[2];
rz(-2.3319862) q[3];
sx q[3];
rz(-1.6428734) q[3];
sx q[3];
rz(-2.0577551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(2.078449) q[0];
sx q[0];
rz(-2.0922631) q[0];
sx q[0];
rz(0.11058841) q[0];
rz(0.94509205) q[1];
sx q[1];
rz(-0.4393622) q[1];
sx q[1];
rz(1.5623215) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54540578) q[0];
sx q[0];
rz(-2.2431787) q[0];
sx q[0];
rz(2.5491712) q[0];
rz(-2.7209372) q[2];
sx q[2];
rz(-1.1569945) q[2];
sx q[2];
rz(-0.54383531) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2436199) q[1];
sx q[1];
rz(-1.0884411) q[1];
sx q[1];
rz(-1.4528758) q[1];
rz(-2.3195295) q[3];
sx q[3];
rz(-1.9371607) q[3];
sx q[3];
rz(-2.0802534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9940146) q[2];
sx q[2];
rz(-0.14573228) q[2];
sx q[2];
rz(0.4501403) q[2];
rz(1.133793) q[3];
sx q[3];
rz(-1.4530051) q[3];
sx q[3];
rz(-0.97761893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56871539) q[0];
sx q[0];
rz(-0.94803888) q[0];
sx q[0];
rz(-0.28133389) q[0];
rz(1.7549134) q[1];
sx q[1];
rz(-1.4298871) q[1];
sx q[1];
rz(2.5414355) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48334405) q[0];
sx q[0];
rz(-1.5243825) q[0];
sx q[0];
rz(0.29057403) q[0];
rz(-pi) q[1];
rz(-1.4900896) q[2];
sx q[2];
rz(-2.5332402) q[2];
sx q[2];
rz(-1.0941015) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.197613) q[1];
sx q[1];
rz(-0.67786067) q[1];
sx q[1];
rz(0.53450905) q[1];
rz(-pi) q[2];
rz(-2.8801877) q[3];
sx q[3];
rz(-0.58901826) q[3];
sx q[3];
rz(2.1847069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0967789) q[2];
sx q[2];
rz(-2.009095) q[2];
sx q[2];
rz(1.0463932) q[2];
rz(-2.7438296) q[3];
sx q[3];
rz(-0.64266509) q[3];
sx q[3];
rz(3.0090295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1716877) q[0];
sx q[0];
rz(-3.1145018) q[0];
sx q[0];
rz(-1.0525674) q[0];
rz(-1.9453847) q[1];
sx q[1];
rz(-1.2095249) q[1];
sx q[1];
rz(0.41950163) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9297816) q[0];
sx q[0];
rz(-2.503184) q[0];
sx q[0];
rz(2.6891461) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2264465) q[2];
sx q[2];
rz(-2.1927877) q[2];
sx q[2];
rz(0.43606191) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4070523) q[1];
sx q[1];
rz(-0.95065763) q[1];
sx q[1];
rz(2.6728515) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8617083) q[3];
sx q[3];
rz(-1.0802764) q[3];
sx q[3];
rz(1.7640997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9516051) q[2];
sx q[2];
rz(-1.1185458) q[2];
sx q[2];
rz(-2.0513963) q[2];
rz(0.53127855) q[3];
sx q[3];
rz(-0.71610206) q[3];
sx q[3];
rz(1.6147015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7200318) q[0];
sx q[0];
rz(-1.5694542) q[0];
sx q[0];
rz(-1.74362) q[0];
rz(0.13294237) q[1];
sx q[1];
rz(-1.839919) q[1];
sx q[1];
rz(-0.001210777) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5842297) q[0];
sx q[0];
rz(-1.5869291) q[0];
sx q[0];
rz(2.7610012) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0693477) q[2];
sx q[2];
rz(-2.8184888) q[2];
sx q[2];
rz(-2.1954775) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1493133) q[1];
sx q[1];
rz(-2.1985801) q[1];
sx q[1];
rz(-1.1924442) q[1];
rz(-pi) q[2];
x q[2];
rz(3.115) q[3];
sx q[3];
rz(-1.4211492) q[3];
sx q[3];
rz(-2.1872471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2289537) q[2];
sx q[2];
rz(-1.8043844) q[2];
sx q[2];
rz(2.9310628) q[2];
rz(2.1052965) q[3];
sx q[3];
rz(-0.784289) q[3];
sx q[3];
rz(-1.9265296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9482816) q[0];
sx q[0];
rz(-1.223215) q[0];
sx q[0];
rz(2.7358828) q[0];
rz(1.7871208) q[1];
sx q[1];
rz(-0.82258075) q[1];
sx q[1];
rz(-2.2192661) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.852823) q[0];
sx q[0];
rz(-1.6772224) q[0];
sx q[0];
rz(1.861051) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0386557) q[2];
sx q[2];
rz(-0.90667001) q[2];
sx q[2];
rz(2.8256106) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3971484) q[1];
sx q[1];
rz(-1.4566629) q[1];
sx q[1];
rz(2.1836917) q[1];
x q[2];
rz(-0.55944632) q[3];
sx q[3];
rz(-1.2370951) q[3];
sx q[3];
rz(-0.73742407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7102082) q[2];
sx q[2];
rz(-1.2146543) q[2];
sx q[2];
rz(-2.7396438) q[2];
rz(-1.2881783) q[3];
sx q[3];
rz(-2.2737019) q[3];
sx q[3];
rz(1.2791546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61889082) q[0];
sx q[0];
rz(-2.0382477) q[0];
sx q[0];
rz(0.064362137) q[0];
rz(-2.3035658) q[1];
sx q[1];
rz(-0.82632724) q[1];
sx q[1];
rz(-1.8109969) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15312402) q[0];
sx q[0];
rz(-2.5543384) q[0];
sx q[0];
rz(0.22756094) q[0];
rz(-pi) q[1];
rz(-0.60637252) q[2];
sx q[2];
rz(-1.1506216) q[2];
sx q[2];
rz(0.5768896) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8992675) q[1];
sx q[1];
rz(-1.0501672) q[1];
sx q[1];
rz(1.8249056) q[1];
rz(-pi) q[2];
rz(-2.182992) q[3];
sx q[3];
rz(-1.7122088) q[3];
sx q[3];
rz(0.32999101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7360721) q[2];
sx q[2];
rz(-1.6330999) q[2];
sx q[2];
rz(2.4580477) q[2];
rz(0.73379597) q[3];
sx q[3];
rz(-1.4132696) q[3];
sx q[3];
rz(0.84764135) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34506327) q[0];
sx q[0];
rz(-1.0836443) q[0];
sx q[0];
rz(1.1731359) q[0];
rz(-0.37551156) q[1];
sx q[1];
rz(-0.48986062) q[1];
sx q[1];
rz(-2.1048022) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11060729) q[0];
sx q[0];
rz(-1.5944423) q[0];
sx q[0];
rz(1.551669) q[0];
x q[1];
rz(-1.8470959) q[2];
sx q[2];
rz(-1.7299998) q[2];
sx q[2];
rz(-0.72180787) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6822676) q[1];
sx q[1];
rz(-2.3258491) q[1];
sx q[1];
rz(1.3774648) q[1];
rz(-pi) q[2];
x q[2];
rz(2.100714) q[3];
sx q[3];
rz(-1.138759) q[3];
sx q[3];
rz(-1.0505067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5658687) q[2];
sx q[2];
rz(-2.2278991) q[2];
sx q[2];
rz(-2.4110528) q[2];
rz(0.20914397) q[3];
sx q[3];
rz(-2.6181965) q[3];
sx q[3];
rz(-1.2701579) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65649477) q[0];
sx q[0];
rz(-0.42335835) q[0];
sx q[0];
rz(0.04743162) q[0];
rz(-0.221953) q[1];
sx q[1];
rz(-2.0675979) q[1];
sx q[1];
rz(2.2304631) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61924261) q[0];
sx q[0];
rz(-1.5372835) q[0];
sx q[0];
rz(1.5438118) q[0];
rz(0.014737562) q[2];
sx q[2];
rz(-0.92916766) q[2];
sx q[2];
rz(0.82889885) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4621689) q[1];
sx q[1];
rz(-2.7738681) q[1];
sx q[1];
rz(-3.1064139) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8857314) q[3];
sx q[3];
rz(-2.5521827) q[3];
sx q[3];
rz(2.3737597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2719416) q[2];
sx q[2];
rz(-2.7615774) q[2];
sx q[2];
rz(-0.16858777) q[2];
rz(0.21909675) q[3];
sx q[3];
rz(-0.85182652) q[3];
sx q[3];
rz(1.8144089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.99437) q[0];
sx q[0];
rz(-0.34039012) q[0];
sx q[0];
rz(2.6841573) q[0];
rz(-1.9153197) q[1];
sx q[1];
rz(-2.8044082) q[1];
sx q[1];
rz(1.7785243) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1049688) q[0];
sx q[0];
rz(-1.5104806) q[0];
sx q[0];
rz(-0.6078267) q[0];
x q[1];
rz(2.8122452) q[2];
sx q[2];
rz(-0.65635704) q[2];
sx q[2];
rz(0.48329566) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7292028) q[1];
sx q[1];
rz(-2.7250184) q[1];
sx q[1];
rz(-0.42804407) q[1];
rz(-pi) q[2];
rz(1.6750338) q[3];
sx q[3];
rz(-2.4564033) q[3];
sx q[3];
rz(2.4327421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7676131) q[2];
sx q[2];
rz(-2.5869936) q[2];
sx q[2];
rz(-2.6895831) q[2];
rz(0.40397817) q[3];
sx q[3];
rz(-1.890506) q[3];
sx q[3];
rz(-2.0154791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6954738) q[0];
sx q[0];
rz(-1.5344545) q[0];
sx q[0];
rz(-2.9547966) q[0];
rz(-2.6192464) q[1];
sx q[1];
rz(-2.6112687) q[1];
sx q[1];
rz(-1.912259) q[1];
rz(-0.28235565) q[2];
sx q[2];
rz(-1.8179802) q[2];
sx q[2];
rz(2.5674934) q[2];
rz(2.8780739) q[3];
sx q[3];
rz(-1.7150039) q[3];
sx q[3];
rz(1.1506469) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
