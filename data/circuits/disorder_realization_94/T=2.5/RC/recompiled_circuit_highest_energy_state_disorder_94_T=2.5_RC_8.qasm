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
rz(-1.0903519) q[0];
sx q[0];
rz(-0.89972377) q[0];
sx q[0];
rz(-1.1296912) q[0];
rz(0.55454412) q[1];
sx q[1];
rz(-0.50538844) q[1];
sx q[1];
rz(-0.55355054) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3239297) q[0];
sx q[0];
rz(-0.50369064) q[0];
sx q[0];
rz(1.5970741) q[0];
rz(-1.1607909) q[2];
sx q[2];
rz(-1.5616817) q[2];
sx q[2];
rz(1.8339008) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.59959902) q[1];
sx q[1];
rz(-2.2570199) q[1];
sx q[1];
rz(-1.047713) q[1];
rz(1.7617202) q[3];
sx q[3];
rz(-0.22293416) q[3];
sx q[3];
rz(0.01252099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0022137) q[2];
sx q[2];
rz(-1.4993818) q[2];
sx q[2];
rz(1.1688983) q[2];
rz(-0.72541952) q[3];
sx q[3];
rz(-1.231671) q[3];
sx q[3];
rz(-1.5554844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91189522) q[0];
sx q[0];
rz(-2.732369) q[0];
sx q[0];
rz(1.0503861) q[0];
rz(-1.3788252) q[1];
sx q[1];
rz(-1.7345112) q[1];
sx q[1];
rz(2.1398267) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79369583) q[0];
sx q[0];
rz(-0.2783723) q[0];
sx q[0];
rz(1.1868366) q[0];
x q[1];
rz(-1.0111079) q[2];
sx q[2];
rz(-0.53072768) q[2];
sx q[2];
rz(-0.84625927) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.92331386) q[1];
sx q[1];
rz(-1.9988195) q[1];
sx q[1];
rz(-1.3175497) q[1];
rz(2.551591) q[3];
sx q[3];
rz(-0.49656061) q[3];
sx q[3];
rz(-2.4499224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.072814552) q[2];
sx q[2];
rz(-0.38663703) q[2];
sx q[2];
rz(-1.6437423) q[2];
rz(-0.98613286) q[3];
sx q[3];
rz(-0.96753263) q[3];
sx q[3];
rz(-0.51924527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3471974) q[0];
sx q[0];
rz(-0.47413844) q[0];
sx q[0];
rz(-0.073632181) q[0];
rz(0.73078784) q[1];
sx q[1];
rz(-1.7981139) q[1];
sx q[1];
rz(2.1582019) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73992621) q[0];
sx q[0];
rz(-2.8186975) q[0];
sx q[0];
rz(-1.0642723) q[0];
rz(-pi) q[1];
rz(2.1165761) q[2];
sx q[2];
rz(-2.3580129) q[2];
sx q[2];
rz(0.5018946) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3427998) q[1];
sx q[1];
rz(-1.4451318) q[1];
sx q[1];
rz(1.8428414) q[1];
x q[2];
rz(-0.077988581) q[3];
sx q[3];
rz(-1.0691563) q[3];
sx q[3];
rz(-2.7671008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.44707766) q[2];
sx q[2];
rz(-0.49459163) q[2];
sx q[2];
rz(3.0818603) q[2];
rz(2.6308502) q[3];
sx q[3];
rz(-0.35664883) q[3];
sx q[3];
rz(-1.3753148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.315595) q[0];
sx q[0];
rz(-0.70846486) q[0];
sx q[0];
rz(-2.3009756) q[0];
rz(2.9261342) q[1];
sx q[1];
rz(-2.0550199) q[1];
sx q[1];
rz(3.0256909) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0784779) q[0];
sx q[0];
rz(-1.3097094) q[0];
sx q[0];
rz(1.9223708) q[0];
rz(-2.5253549) q[2];
sx q[2];
rz(-0.35348693) q[2];
sx q[2];
rz(-0.36187672) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2825968) q[1];
sx q[1];
rz(-0.94817417) q[1];
sx q[1];
rz(2.7864561) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41036116) q[3];
sx q[3];
rz(-1.6660353) q[3];
sx q[3];
rz(1.8382753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0243715) q[2];
sx q[2];
rz(-1.8481959) q[2];
sx q[2];
rz(2.2852066) q[2];
rz(-2.3971426) q[3];
sx q[3];
rz(-0.88862935) q[3];
sx q[3];
rz(-2.7763247) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15904467) q[0];
sx q[0];
rz(-0.58369517) q[0];
sx q[0];
rz(-0.99345508) q[0];
rz(-2.0221201) q[1];
sx q[1];
rz(-1.9529724) q[1];
sx q[1];
rz(3.1192034) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7763045) q[0];
sx q[0];
rz(-2.9927706) q[0];
sx q[0];
rz(-2.5689199) q[0];
x q[1];
rz(-2.585942) q[2];
sx q[2];
rz(-1.0684864) q[2];
sx q[2];
rz(2.7567418) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7931079) q[1];
sx q[1];
rz(-0.99103084) q[1];
sx q[1];
rz(-2.0391885) q[1];
x q[2];
rz(-2.5804912) q[3];
sx q[3];
rz(-2.3659539) q[3];
sx q[3];
rz(-0.55828324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8523954) q[2];
sx q[2];
rz(-0.22326938) q[2];
sx q[2];
rz(-1.7247464) q[2];
rz(1.019574) q[3];
sx q[3];
rz(-2.151078) q[3];
sx q[3];
rz(2.8250601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.9483017) q[0];
sx q[0];
rz(-3.1315849) q[0];
sx q[0];
rz(-0.77578068) q[0];
rz(-1.7316445) q[1];
sx q[1];
rz(-1.3103139) q[1];
sx q[1];
rz(1.5595248) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55627953) q[0];
sx q[0];
rz(-1.5405498) q[0];
sx q[0];
rz(-1.6325238) q[0];
x q[1];
rz(2.7804017) q[2];
sx q[2];
rz(-1.0345999) q[2];
sx q[2];
rz(2.9597832) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3969793) q[1];
sx q[1];
rz(-1.1049264) q[1];
sx q[1];
rz(-0.14772285) q[1];
x q[2];
rz(-1.9385064) q[3];
sx q[3];
rz(-1.0818521) q[3];
sx q[3];
rz(0.5209825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5887624) q[2];
sx q[2];
rz(-1.7071807) q[2];
sx q[2];
rz(-0.32858783) q[2];
rz(1.5594679) q[3];
sx q[3];
rz(-0.72454536) q[3];
sx q[3];
rz(2.4174378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.5174823) q[0];
sx q[0];
rz(-0.8172074) q[0];
sx q[0];
rz(-0.66993129) q[0];
rz(1.6239369) q[1];
sx q[1];
rz(-2.0567963) q[1];
sx q[1];
rz(-2.0673015) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6321323) q[0];
sx q[0];
rz(-1.355251) q[0];
sx q[0];
rz(2.4556405) q[0];
x q[1];
rz(-0.36342095) q[2];
sx q[2];
rz(-1.7426024) q[2];
sx q[2];
rz(-1.8653143) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5902293) q[1];
sx q[1];
rz(-1.6632212) q[1];
sx q[1];
rz(-1.2552099) q[1];
x q[2];
rz(-2.9250894) q[3];
sx q[3];
rz(-2.2364383) q[3];
sx q[3];
rz(0.82279278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.39277789) q[2];
sx q[2];
rz(-0.98069507) q[2];
sx q[2];
rz(-1.3856212) q[2];
rz(-0.70220002) q[3];
sx q[3];
rz(-2.8148013) q[3];
sx q[3];
rz(2.0328111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7053407) q[0];
sx q[0];
rz(-1.5441283) q[0];
sx q[0];
rz(-2.3633549) q[0];
rz(-1.2938007) q[1];
sx q[1];
rz(-1.8840645) q[1];
sx q[1];
rz(0.055880849) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.125408) q[0];
sx q[0];
rz(-1.8620389) q[0];
sx q[0];
rz(1.6470558) q[0];
x q[1];
rz(0.22578366) q[2];
sx q[2];
rz(-1.7541382) q[2];
sx q[2];
rz(2.8933337) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3371007) q[1];
sx q[1];
rz(-2.7332515) q[1];
sx q[1];
rz(-0.33894811) q[1];
rz(-2.6229739) q[3];
sx q[3];
rz(-2.2317776) q[3];
sx q[3];
rz(0.80693124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.087626584) q[2];
sx q[2];
rz(-1.0216594) q[2];
sx q[2];
rz(-0.6915687) q[2];
rz(-0.21928731) q[3];
sx q[3];
rz(-1.6577474) q[3];
sx q[3];
rz(-1.9482435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9241065) q[0];
sx q[0];
rz(-3.1052488) q[0];
sx q[0];
rz(-2.0935667) q[0];
rz(-0.74608392) q[1];
sx q[1];
rz(-1.289295) q[1];
sx q[1];
rz(-2.7229436) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0264473) q[0];
sx q[0];
rz(-1.5799684) q[0];
sx q[0];
rz(-1.5546726) q[0];
rz(-0.92709686) q[2];
sx q[2];
rz(-1.3852556) q[2];
sx q[2];
rz(0.58840051) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9050952) q[1];
sx q[1];
rz(-1.4346544) q[1];
sx q[1];
rz(0.13827583) q[1];
rz(-pi) q[2];
rz(-0.04808406) q[3];
sx q[3];
rz(-0.86838956) q[3];
sx q[3];
rz(2.4772205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6333542) q[2];
sx q[2];
rz(-2.5373122) q[2];
sx q[2];
rz(-2.5227127) q[2];
rz(-0.45281705) q[3];
sx q[3];
rz(-1.6524977) q[3];
sx q[3];
rz(0.49310163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5433237) q[0];
sx q[0];
rz(-3.0040574) q[0];
sx q[0];
rz(-1.8813053) q[0];
rz(-2.8434985) q[1];
sx q[1];
rz(-2.1838078) q[1];
sx q[1];
rz(1.4990998) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84628203) q[0];
sx q[0];
rz(-0.75937004) q[0];
sx q[0];
rz(1.6037386) q[0];
rz(-pi) q[1];
rz(1.4269198) q[2];
sx q[2];
rz(-1.292079) q[2];
sx q[2];
rz(0.65482322) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.19263515) q[1];
sx q[1];
rz(-1.7970553) q[1];
sx q[1];
rz(0.016955777) q[1];
x q[2];
rz(2.9931287) q[3];
sx q[3];
rz(-1.1824596) q[3];
sx q[3];
rz(-2.9290105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.62557158) q[2];
sx q[2];
rz(-1.2837774) q[2];
sx q[2];
rz(0.18216356) q[2];
rz(-2.9158909) q[3];
sx q[3];
rz(-2.1322026) q[3];
sx q[3];
rz(1.766248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43212776) q[0];
sx q[0];
rz(-1.3668677) q[0];
sx q[0];
rz(0.9912542) q[0];
rz(-1.9215012) q[1];
sx q[1];
rz(-2.5972001) q[1];
sx q[1];
rz(-0.77144815) q[1];
rz(2.0408199) q[2];
sx q[2];
rz(-0.19954544) q[2];
sx q[2];
rz(-0.96878845) q[2];
rz(2.3516918) q[3];
sx q[3];
rz(-0.98179437) q[3];
sx q[3];
rz(-2.475987) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
