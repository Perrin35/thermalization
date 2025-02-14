OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9187501) q[0];
sx q[0];
rz(-1.1981244) q[0];
sx q[0];
rz(-2.591748) q[0];
rz(-0.066601872) q[1];
sx q[1];
rz(5.0636518) q[1];
sx q[1];
rz(10.650462) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8526575) q[0];
sx q[0];
rz(-2.9338917) q[0];
sx q[0];
rz(3.0124979) q[0];
x q[1];
rz(0.066402175) q[2];
sx q[2];
rz(-1.7207256) q[2];
sx q[2];
rz(-0.87519803) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1154626) q[1];
sx q[1];
rz(-1.5688174) q[1];
sx q[1];
rz(1.8972562) q[1];
x q[2];
rz(0.59210316) q[3];
sx q[3];
rz(-1.0463011) q[3];
sx q[3];
rz(-1.2712353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.77039069) q[2];
sx q[2];
rz(-1.0545701) q[2];
sx q[2];
rz(-2.2626109) q[2];
rz(-0.061035872) q[3];
sx q[3];
rz(-1.4417442) q[3];
sx q[3];
rz(2.216831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9691129) q[0];
sx q[0];
rz(-1.1107439) q[0];
sx q[0];
rz(-1.9222395) q[0];
rz(-2.1439233) q[1];
sx q[1];
rz(-1.6976633) q[1];
sx q[1];
rz(0.77233058) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15914842) q[0];
sx q[0];
rz(-1.8799025) q[0];
sx q[0];
rz(-2.3594666) q[0];
rz(1.6222811) q[2];
sx q[2];
rz(-0.33308187) q[2];
sx q[2];
rz(-0.55120984) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0454084) q[1];
sx q[1];
rz(-1.533877) q[1];
sx q[1];
rz(0.30942076) q[1];
rz(-pi) q[2];
rz(3.1082013) q[3];
sx q[3];
rz(-0.62761939) q[3];
sx q[3];
rz(-0.96783584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4274009) q[2];
sx q[2];
rz(-1.4560459) q[2];
sx q[2];
rz(-1.5217155) q[2];
rz(1.4766258) q[3];
sx q[3];
rz(-1.0790389) q[3];
sx q[3];
rz(2.9286706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.49887) q[0];
sx q[0];
rz(-1.9864137) q[0];
sx q[0];
rz(2.2464519) q[0];
rz(0.25998947) q[1];
sx q[1];
rz(-1.7951868) q[1];
sx q[1];
rz(2.3666429) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6942368) q[0];
sx q[0];
rz(-1.0123582) q[0];
sx q[0];
rz(1.9959227) q[0];
x q[1];
rz(2.9812212) q[2];
sx q[2];
rz(-0.43284349) q[2];
sx q[2];
rz(0.835604) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4878825) q[1];
sx q[1];
rz(-2.0984432) q[1];
sx q[1];
rz(0.58077537) q[1];
rz(-pi) q[2];
rz(0.633274) q[3];
sx q[3];
rz(-0.97304854) q[3];
sx q[3];
rz(0.55766314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.2669175) q[2];
sx q[2];
rz(-2.3053034) q[2];
sx q[2];
rz(2.3036892) q[2];
rz(-0.72495929) q[3];
sx q[3];
rz(-0.91491142) q[3];
sx q[3];
rz(-0.27423492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2570141) q[0];
sx q[0];
rz(-3.0138636) q[0];
sx q[0];
rz(1.9789486) q[0];
rz(-1.1391901) q[1];
sx q[1];
rz(-1.2290686) q[1];
sx q[1];
rz(0.3831648) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98477798) q[0];
sx q[0];
rz(-1.1582644) q[0];
sx q[0];
rz(-1.5810313) q[0];
rz(1.2960984) q[2];
sx q[2];
rz(-0.75680671) q[2];
sx q[2];
rz(-1.5001378) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5096915) q[1];
sx q[1];
rz(-1.5576524) q[1];
sx q[1];
rz(0.34850328) q[1];
rz(-pi) q[2];
rz(1.9145033) q[3];
sx q[3];
rz(-2.1919804) q[3];
sx q[3];
rz(-0.4200926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7202516) q[2];
sx q[2];
rz(-0.64695078) q[2];
sx q[2];
rz(1.425239) q[2];
rz(-2.0187812) q[3];
sx q[3];
rz(-2.4425826) q[3];
sx q[3];
rz(0.56306806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7177893) q[0];
sx q[0];
rz(-0.20861067) q[0];
sx q[0];
rz(-2.8243689) q[0];
rz(-0.18103655) q[1];
sx q[1];
rz(-1.92417) q[1];
sx q[1];
rz(-2.5203629) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3040627) q[0];
sx q[0];
rz(-0.777839) q[0];
sx q[0];
rz(-1.8767979) q[0];
rz(-pi) q[1];
x q[1];
rz(0.076940342) q[2];
sx q[2];
rz(-2.6443993) q[2];
sx q[2];
rz(0.92958462) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9108429) q[1];
sx q[1];
rz(-0.98346868) q[1];
sx q[1];
rz(3.0452646) q[1];
rz(-pi) q[2];
rz(2.0192573) q[3];
sx q[3];
rz(-2.3178007) q[3];
sx q[3];
rz(-0.5169249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.3250331) q[2];
sx q[2];
rz(-0.72924048) q[2];
sx q[2];
rz(1.0664335) q[2];
rz(-2.1224497) q[3];
sx q[3];
rz(-2.416553) q[3];
sx q[3];
rz(-1.2044005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8891193) q[0];
sx q[0];
rz(-0.42943615) q[0];
sx q[0];
rz(0.47527894) q[0];
rz(-2.1976082) q[1];
sx q[1];
rz(-1.9119268) q[1];
sx q[1];
rz(-0.42517391) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89579534) q[0];
sx q[0];
rz(-2.6962792) q[0];
sx q[0];
rz(-2.6209339) q[0];
rz(-pi) q[1];
rz(2.6475372) q[2];
sx q[2];
rz(-1.5580342) q[2];
sx q[2];
rz(-0.27496613) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2819967) q[1];
sx q[1];
rz(-0.55742555) q[1];
sx q[1];
rz(-1.2453399) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75192368) q[3];
sx q[3];
rz(-2.6787191) q[3];
sx q[3];
rz(1.9583256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7298594) q[2];
sx q[2];
rz(-1.6526165) q[2];
sx q[2];
rz(-2.497351) q[2];
rz(-1.980137) q[3];
sx q[3];
rz(-0.96122733) q[3];
sx q[3];
rz(0.78288356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92530584) q[0];
sx q[0];
rz(-1.1290978) q[0];
sx q[0];
rz(0.5710477) q[0];
rz(1.1043999) q[1];
sx q[1];
rz(-2.0490502) q[1];
sx q[1];
rz(0.33448321) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0153941) q[0];
sx q[0];
rz(-1.6078464) q[0];
sx q[0];
rz(-1.1657715) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38774298) q[2];
sx q[2];
rz(-0.84008145) q[2];
sx q[2];
rz(0.71114388) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.7664288) q[1];
sx q[1];
rz(-0.42459449) q[1];
sx q[1];
rz(-2.2284075) q[1];
rz(-pi) q[2];
rz(-1.2932106) q[3];
sx q[3];
rz(-2.2065139) q[3];
sx q[3];
rz(-1.1073974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.107848) q[2];
sx q[2];
rz(-1.1128384) q[2];
sx q[2];
rz(0.10438485) q[2];
rz(0.32621041) q[3];
sx q[3];
rz(-0.86745894) q[3];
sx q[3];
rz(2.4645658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0584843) q[0];
sx q[0];
rz(-1.2106189) q[0];
sx q[0];
rz(0.50115681) q[0];
rz(-1.1309364) q[1];
sx q[1];
rz(-0.94291818) q[1];
sx q[1];
rz(-1.8556192) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0752397) q[0];
sx q[0];
rz(-1.4202002) q[0];
sx q[0];
rz(2.0429918) q[0];
rz(-pi) q[1];
rz(2.7825955) q[2];
sx q[2];
rz(-2.411826) q[2];
sx q[2];
rz(-1.9791918) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0803538) q[1];
sx q[1];
rz(-0.29268943) q[1];
sx q[1];
rz(0.55673843) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2896832) q[3];
sx q[3];
rz(-2.503667) q[3];
sx q[3];
rz(-1.7450116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8274294) q[2];
sx q[2];
rz(-2.8848727) q[2];
sx q[2];
rz(2.2641505) q[2];
rz(-0.99831239) q[3];
sx q[3];
rz(-0.58019296) q[3];
sx q[3];
rz(-1.4210526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(1.5494635) q[0];
sx q[0];
rz(-1.9100459) q[0];
sx q[0];
rz(-0.25729427) q[0];
rz(-0.65722242) q[1];
sx q[1];
rz(-2.6723599) q[1];
sx q[1];
rz(1.8990272) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3982733) q[0];
sx q[0];
rz(-0.87617517) q[0];
sx q[0];
rz(2.1477703) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3506566) q[2];
sx q[2];
rz(-1.4688014) q[2];
sx q[2];
rz(-1.9079218) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4081289) q[1];
sx q[1];
rz(-2.4545547) q[1];
sx q[1];
rz(-2.4279159) q[1];
x q[2];
rz(2.8065422) q[3];
sx q[3];
rz(-1.3962702) q[3];
sx q[3];
rz(-1.9843335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.88951748) q[2];
sx q[2];
rz(-1.4318848) q[2];
sx q[2];
rz(2.4375708) q[2];
rz(1.6440803) q[3];
sx q[3];
rz(-1.5181395) q[3];
sx q[3];
rz(1.2161072) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4154476) q[0];
sx q[0];
rz(-2.4496267) q[0];
sx q[0];
rz(-2.6494001) q[0];
rz(-0.65912143) q[1];
sx q[1];
rz(-0.4159795) q[1];
sx q[1];
rz(-0.98512828) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8644442) q[0];
sx q[0];
rz(-1.6246968) q[0];
sx q[0];
rz(3.0646851) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4617301) q[2];
sx q[2];
rz(-1.1972053) q[2];
sx q[2];
rz(1.2078169) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.0041005153) q[1];
sx q[1];
rz(-1.8086063) q[1];
sx q[1];
rz(0.40356839) q[1];
x q[2];
rz(2.3045818) q[3];
sx q[3];
rz(-1.600255) q[3];
sx q[3];
rz(-3.000976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8681086) q[2];
sx q[2];
rz(-2.2874139) q[2];
sx q[2];
rz(-1.9197397) q[2];
rz(-2.6814804) q[3];
sx q[3];
rz(-1.8729788) q[3];
sx q[3];
rz(-2.4251895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39501326) q[0];
sx q[0];
rz(-1.6759251) q[0];
sx q[0];
rz(1.2431086) q[0];
rz(-1.5383491) q[1];
sx q[1];
rz(-2.1080882) q[1];
sx q[1];
rz(2.7079667) q[1];
rz(-1.4547841) q[2];
sx q[2];
rz(-2.3064936) q[2];
sx q[2];
rz(2.4207122) q[2];
rz(-0.015751377) q[3];
sx q[3];
rz(-1.8498265) q[3];
sx q[3];
rz(-2.7791666) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
