OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71830463) q[0];
sx q[0];
rz(-0.77646065) q[0];
sx q[0];
rz(2.7756696) q[0];
rz(-2.794682) q[1];
sx q[1];
rz(3.4263098) q[1];
sx q[1];
rz(10.813536) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4735218) q[0];
sx q[0];
rz(-1.6374303) q[0];
sx q[0];
rz(-0.20247831) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69268815) q[2];
sx q[2];
rz(-1.8028304) q[2];
sx q[2];
rz(1.2218066) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5386913) q[1];
sx q[1];
rz(-0.60990342) q[1];
sx q[1];
rz(-2.2944488) q[1];
rz(-0.71309375) q[3];
sx q[3];
rz(-0.72848195) q[3];
sx q[3];
rz(0.99458867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.95691386) q[2];
sx q[2];
rz(-2.2214486) q[2];
sx q[2];
rz(-0.93519768) q[2];
rz(-2.8397371) q[3];
sx q[3];
rz(-1.5993092) q[3];
sx q[3];
rz(2.7479318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.1644156) q[0];
sx q[0];
rz(-1.2657607) q[0];
sx q[0];
rz(-2.6432977) q[0];
rz(-1.4277108) q[1];
sx q[1];
rz(-1.9352813) q[1];
sx q[1];
rz(2.0548342) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55288494) q[0];
sx q[0];
rz(-2.4037139) q[0];
sx q[0];
rz(0.35564519) q[0];
rz(-pi) q[1];
rz(-1.3042167) q[2];
sx q[2];
rz(-1.3836897) q[2];
sx q[2];
rz(3.0361036) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1732498) q[1];
sx q[1];
rz(-1.1361215) q[1];
sx q[1];
rz(1.7477075) q[1];
rz(-pi) q[2];
rz(1.5828176) q[3];
sx q[3];
rz(-1.2526027) q[3];
sx q[3];
rz(-2.1592922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9429417) q[2];
sx q[2];
rz(-2.8309839) q[2];
sx q[2];
rz(-1.8093713) q[2];
rz(1.9940935) q[3];
sx q[3];
rz(-1.56286) q[3];
sx q[3];
rz(-1.3326299) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20951095) q[0];
sx q[0];
rz(-1.4028343) q[0];
sx q[0];
rz(-2.984356) q[0];
rz(-1.2278185) q[1];
sx q[1];
rz(-0.20918748) q[1];
sx q[1];
rz(-0.42207178) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7868766) q[0];
sx q[0];
rz(-1.3930969) q[0];
sx q[0];
rz(0.091376809) q[0];
x q[1];
rz(-2.5890805) q[2];
sx q[2];
rz(-0.32304155) q[2];
sx q[2];
rz(-0.99977899) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8542618) q[1];
sx q[1];
rz(-0.39188436) q[1];
sx q[1];
rz(1.803483) q[1];
rz(-1.8881599) q[3];
sx q[3];
rz(-2.3471213) q[3];
sx q[3];
rz(1.0004071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.32855365) q[2];
sx q[2];
rz(-0.96031323) q[2];
sx q[2];
rz(-0.73406827) q[2];
rz(-2.0640533) q[3];
sx q[3];
rz(-2.4534093) q[3];
sx q[3];
rz(3.0663826) q[3];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5587191) q[0];
sx q[0];
rz(-2.131077) q[0];
sx q[0];
rz(-3.1251113) q[0];
rz(3.0670498) q[1];
sx q[1];
rz(-2.6838979) q[1];
sx q[1];
rz(0.25513908) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3203188) q[0];
sx q[0];
rz(-0.50423586) q[0];
sx q[0];
rz(-2.8404499) q[0];
rz(-2.1367777) q[2];
sx q[2];
rz(-1.2850956) q[2];
sx q[2];
rz(0.20488747) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.29683477) q[1];
sx q[1];
rz(-0.95176178) q[1];
sx q[1];
rz(2.1059787) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0276518) q[3];
sx q[3];
rz(-1.9181532) q[3];
sx q[3];
rz(1.243227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.61465803) q[2];
sx q[2];
rz(-0.69710985) q[2];
sx q[2];
rz(-2.441791) q[2];
rz(0.091863306) q[3];
sx q[3];
rz(-1.2173165) q[3];
sx q[3];
rz(0.45472586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5192473) q[0];
sx q[0];
rz(-0.82038251) q[0];
sx q[0];
rz(2.0297594) q[0];
rz(-0.19019292) q[1];
sx q[1];
rz(-2.2564502) q[1];
sx q[1];
rz(-1.9817188) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58901133) q[0];
sx q[0];
rz(-1.9144626) q[0];
sx q[0];
rz(-0.036618311) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14638193) q[2];
sx q[2];
rz(-0.9375776) q[2];
sx q[2];
rz(-0.13356471) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3768757) q[1];
sx q[1];
rz(-2.5333166) q[1];
sx q[1];
rz(2.6981955) q[1];
rz(1.0489063) q[3];
sx q[3];
rz(-1.5318222) q[3];
sx q[3];
rz(-0.029595395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5707034) q[2];
sx q[2];
rz(-2.3184226) q[2];
sx q[2];
rz(2.9111351) q[2];
rz(-1.0620091) q[3];
sx q[3];
rz(-1.8657203) q[3];
sx q[3];
rz(-1.6331858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.688711) q[0];
sx q[0];
rz(-1.0588366) q[0];
sx q[0];
rz(-1.3558615) q[0];
rz(0.52976766) q[1];
sx q[1];
rz(-0.7822839) q[1];
sx q[1];
rz(-1.1518325) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065581948) q[0];
sx q[0];
rz(-1.9261203) q[0];
sx q[0];
rz(1.6760582) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.095345796) q[2];
sx q[2];
rz(-1.0311677) q[2];
sx q[2];
rz(-2.6971872) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.63910145) q[1];
sx q[1];
rz(-1.3768973) q[1];
sx q[1];
rz(-0.37096937) q[1];
rz(-0.21422503) q[3];
sx q[3];
rz(-1.6856632) q[3];
sx q[3];
rz(1.3319931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4749703) q[2];
sx q[2];
rz(-0.64133659) q[2];
sx q[2];
rz(-2.6744911) q[2];
rz(-1.0111672) q[3];
sx q[3];
rz(-0.34365383) q[3];
sx q[3];
rz(-1.3694192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8491299) q[0];
sx q[0];
rz(-0.15501538) q[0];
sx q[0];
rz(-0.22468654) q[0];
rz(-2.977773) q[1];
sx q[1];
rz(-0.33912173) q[1];
sx q[1];
rz(-0.64635578) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77262473) q[0];
sx q[0];
rz(-0.59212055) q[0];
sx q[0];
rz(-1.1117685) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1557344) q[2];
sx q[2];
rz(-0.38160593) q[2];
sx q[2];
rz(-1.3473036) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7228292) q[1];
sx q[1];
rz(-1.495528) q[1];
sx q[1];
rz(-2.597517) q[1];
rz(-pi) q[2];
rz(0.38742806) q[3];
sx q[3];
rz(-1.6968075) q[3];
sx q[3];
rz(2.2962017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9859163) q[2];
sx q[2];
rz(-1.4480269) q[2];
sx q[2];
rz(-0.4099561) q[2];
rz(-2.3112467) q[3];
sx q[3];
rz(-2.1928619) q[3];
sx q[3];
rz(1.4478987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(2.603867) q[0];
sx q[0];
rz(-2.4483838) q[0];
sx q[0];
rz(-0.24630462) q[0];
rz(2.314997) q[1];
sx q[1];
rz(-1.7431424) q[1];
sx q[1];
rz(-2.8776317) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9519099) q[0];
sx q[0];
rz(-3.0458458) q[0];
sx q[0];
rz(-1.4950947) q[0];
rz(0.47672456) q[2];
sx q[2];
rz(-0.84273224) q[2];
sx q[2];
rz(-1.5155033) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3726793) q[1];
sx q[1];
rz(-2.2376334) q[1];
sx q[1];
rz(-2.7808378) q[1];
rz(-pi) q[2];
rz(2.554515) q[3];
sx q[3];
rz(-1.2012606) q[3];
sx q[3];
rz(-2.0943506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.10494122) q[2];
sx q[2];
rz(-1.2678601) q[2];
sx q[2];
rz(-2.4264753) q[2];
rz(-3.0702843) q[3];
sx q[3];
rz(-1.5569867) q[3];
sx q[3];
rz(-2.3947072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0095373) q[0];
sx q[0];
rz(-1.4893463) q[0];
sx q[0];
rz(-0.43553964) q[0];
rz(-0.14777331) q[1];
sx q[1];
rz(-1.4316033) q[1];
sx q[1];
rz(-1.6730283) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4184121) q[0];
sx q[0];
rz(-1.3910798) q[0];
sx q[0];
rz(2.7030858) q[0];
x q[1];
rz(2.4442102) q[2];
sx q[2];
rz(-1.6314403) q[2];
sx q[2];
rz(1.9549119) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.802001) q[1];
sx q[1];
rz(-2.7164384) q[1];
sx q[1];
rz(1.156686) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7870258) q[3];
sx q[3];
rz(-1.9659967) q[3];
sx q[3];
rz(2.59557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.74464166) q[2];
sx q[2];
rz(-1.3775185) q[2];
sx q[2];
rz(-0.89293876) q[2];
rz(1.2785771) q[3];
sx q[3];
rz(-1.1568926) q[3];
sx q[3];
rz(-1.4634092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-2.8466012) q[0];
sx q[0];
rz(-2.8104267) q[0];
sx q[0];
rz(-0.60605979) q[0];
rz(3.1312969) q[1];
sx q[1];
rz(-0.98774397) q[1];
sx q[1];
rz(-3.0616679) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29045668) q[0];
sx q[0];
rz(-0.62617597) q[0];
sx q[0];
rz(2.0100223) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8907919) q[2];
sx q[2];
rz(-0.75437949) q[2];
sx q[2];
rz(0.88448521) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8070302) q[1];
sx q[1];
rz(-1.9777781) q[1];
sx q[1];
rz(1.5949859) q[1];
rz(-0.20584835) q[3];
sx q[3];
rz(-1.7966812) q[3];
sx q[3];
rz(-0.0061279002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.14110485) q[2];
sx q[2];
rz(-2.2147369) q[2];
sx q[2];
rz(-1.0151939) q[2];
rz(-0.60123932) q[3];
sx q[3];
rz(-2.4398118) q[3];
sx q[3];
rz(-2.0950441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4973608) q[0];
sx q[0];
rz(-1.5174706) q[0];
sx q[0];
rz(0.68751412) q[0];
rz(2.5524706) q[1];
sx q[1];
rz(-2.4644869) q[1];
sx q[1];
rz(-0.45225515) q[1];
rz(-1.4745787) q[2];
sx q[2];
rz(-2.0141891) q[2];
sx q[2];
rz(-0.63465848) q[2];
rz(0.11832931) q[3];
sx q[3];
rz(-0.82810267) q[3];
sx q[3];
rz(-1.8406263) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
