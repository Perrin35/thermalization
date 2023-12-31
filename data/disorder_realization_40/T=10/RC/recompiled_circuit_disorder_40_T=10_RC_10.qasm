OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6089132) q[0];
sx q[0];
rz(-0.37663868) q[0];
sx q[0];
rz(-3.0298046) q[0];
rz(1.6821661) q[1];
sx q[1];
rz(-1.4844126) q[1];
sx q[1];
rz(2.9878374) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35225866) q[0];
sx q[0];
rz(-2.6129122) q[0];
sx q[0];
rz(-2.5369011) q[0];
rz(-pi) q[1];
x q[1];
rz(0.076924952) q[2];
sx q[2];
rz(-1.3300606) q[2];
sx q[2];
rz(-1.392729) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.54290463) q[1];
sx q[1];
rz(-0.3202657) q[1];
sx q[1];
rz(-1.4742875) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1313685) q[3];
sx q[3];
rz(-1.2636375) q[3];
sx q[3];
rz(2.9205521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.443632) q[2];
sx q[2];
rz(-1.4322832) q[2];
sx q[2];
rz(1.4367746) q[2];
rz(-0.73389655) q[3];
sx q[3];
rz(-1.5926444) q[3];
sx q[3];
rz(0.51600391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.88965082) q[0];
sx q[0];
rz(-1.9152859) q[0];
sx q[0];
rz(2.2170128) q[0];
rz(2.1444767) q[1];
sx q[1];
rz(-2.6328502) q[1];
sx q[1];
rz(1.3234214) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91081496) q[0];
sx q[0];
rz(-1.9998904) q[0];
sx q[0];
rz(0.50165117) q[0];
rz(-pi) q[1];
rz(-0.3496062) q[2];
sx q[2];
rz(-1.6125624) q[2];
sx q[2];
rz(1.8984399) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1393226) q[1];
sx q[1];
rz(-1.8404669) q[1];
sx q[1];
rz(-2.3766999) q[1];
x q[2];
rz(-1.0074535) q[3];
sx q[3];
rz(-0.65348071) q[3];
sx q[3];
rz(0.37936488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.20415846) q[2];
sx q[2];
rz(-1.5896475) q[2];
sx q[2];
rz(0.75817529) q[2];
rz(-2.5126863) q[3];
sx q[3];
rz(-2.7401676) q[3];
sx q[3];
rz(-1.988407) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73308289) q[0];
sx q[0];
rz(-1.2671616) q[0];
sx q[0];
rz(0.68840233) q[0];
rz(-3.0738661) q[1];
sx q[1];
rz(-1.7522782) q[1];
sx q[1];
rz(0.53007954) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2358658) q[0];
sx q[0];
rz(-1.6775963) q[0];
sx q[0];
rz(-1.8624767) q[0];
x q[1];
rz(2.2377551) q[2];
sx q[2];
rz(-2.505216) q[2];
sx q[2];
rz(2.9449376) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.053777545) q[1];
sx q[1];
rz(-2.3936845) q[1];
sx q[1];
rz(-2.0762216) q[1];
rz(-pi) q[2];
rz(-3.043626) q[3];
sx q[3];
rz(-1.8091822) q[3];
sx q[3];
rz(2.1408248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.80785859) q[2];
sx q[2];
rz(-0.01161751) q[2];
sx q[2];
rz(-0.90144908) q[2];
rz(-2.3060913) q[3];
sx q[3];
rz(-1.52799) q[3];
sx q[3];
rz(1.8301331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9505342) q[0];
sx q[0];
rz(-1.598851) q[0];
sx q[0];
rz(-0.78432551) q[0];
rz(-0.061231881) q[1];
sx q[1];
rz(-2.4274554) q[1];
sx q[1];
rz(-3.004946) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.559583) q[0];
sx q[0];
rz(-2.0499381) q[0];
sx q[0];
rz(-0.14736202) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58089528) q[2];
sx q[2];
rz(-0.53933203) q[2];
sx q[2];
rz(2.5748594) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1130484) q[1];
sx q[1];
rz(-1.5469157) q[1];
sx q[1];
rz(-1.5024115) q[1];
x q[2];
rz(0.99478787) q[3];
sx q[3];
rz(-1.4928865) q[3];
sx q[3];
rz(-0.60914492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0121997) q[2];
sx q[2];
rz(-0.94038525) q[2];
sx q[2];
rz(0.56048918) q[2];
rz(0.012332049) q[3];
sx q[3];
rz(-0.90356946) q[3];
sx q[3];
rz(2.0509317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43301582) q[0];
sx q[0];
rz(-2.5362159) q[0];
sx q[0];
rz(0.82114712) q[0];
rz(-2.2654146) q[1];
sx q[1];
rz(-0.89996243) q[1];
sx q[1];
rz(1.7339773) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27027425) q[0];
sx q[0];
rz(-0.5373913) q[0];
sx q[0];
rz(2.4094765) q[0];
rz(2.0166964) q[2];
sx q[2];
rz(-0.6302399) q[2];
sx q[2];
rz(-1.557204) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.447532) q[1];
sx q[1];
rz(-0.58528712) q[1];
sx q[1];
rz(-1.3259757) q[1];
rz(-pi) q[2];
rz(1.9305265) q[3];
sx q[3];
rz(-1.3621646) q[3];
sx q[3];
rz(1.2071351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.451482) q[2];
sx q[2];
rz(-1.2146981) q[2];
sx q[2];
rz(-0.042479854) q[2];
rz(2.5111607) q[3];
sx q[3];
rz(-0.63215956) q[3];
sx q[3];
rz(0.49155864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38934389) q[0];
sx q[0];
rz(-1.9135973) q[0];
sx q[0];
rz(-2.65843) q[0];
rz(1.0522316) q[1];
sx q[1];
rz(-1.9960884) q[1];
sx q[1];
rz(2.5767456) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72950596) q[0];
sx q[0];
rz(-0.54486638) q[0];
sx q[0];
rz(1.3077523) q[0];
rz(1.9916612) q[2];
sx q[2];
rz(-1.5016342) q[2];
sx q[2];
rz(1.3644621) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40183345) q[1];
sx q[1];
rz(-1.9976915) q[1];
sx q[1];
rz(-0.065211936) q[1];
x q[2];
rz(-1.4156028) q[3];
sx q[3];
rz(-1.3421913) q[3];
sx q[3];
rz(2.5582563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.57006449) q[2];
sx q[2];
rz(-2.0740985) q[2];
sx q[2];
rz(1.338039) q[2];
rz(-1.8367052) q[3];
sx q[3];
rz(-2.0740502) q[3];
sx q[3];
rz(-0.6750955) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17334443) q[0];
sx q[0];
rz(-1.4302379) q[0];
sx q[0];
rz(2.5937953) q[0];
rz(0.785218) q[1];
sx q[1];
rz(-1.3350057) q[1];
sx q[1];
rz(2.8731667) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1584394) q[0];
sx q[0];
rz(-1.325325) q[0];
sx q[0];
rz(-1.655274) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3556049) q[2];
sx q[2];
rz(-2.18835) q[2];
sx q[2];
rz(0.76030375) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5749579) q[1];
sx q[1];
rz(-2.1556427) q[1];
sx q[1];
rz(-3.023874) q[1];
x q[2];
rz(-1.0309585) q[3];
sx q[3];
rz(-1.8317878) q[3];
sx q[3];
rz(0.48418448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5238374) q[2];
sx q[2];
rz(-0.80344168) q[2];
sx q[2];
rz(-0.8141554) q[2];
rz(2.7653149) q[3];
sx q[3];
rz(-1.9777931) q[3];
sx q[3];
rz(-3.0686839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5557264) q[0];
sx q[0];
rz(-1.7365475) q[0];
sx q[0];
rz(-1.0193753) q[0];
rz(-2.2881919) q[1];
sx q[1];
rz(-1.1420206) q[1];
sx q[1];
rz(-0.44874915) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42496029) q[0];
sx q[0];
rz(-0.53629959) q[0];
sx q[0];
rz(-2.0400356) q[0];
rz(-pi) q[1];
rz(0.25787392) q[2];
sx q[2];
rz(-1.2784625) q[2];
sx q[2];
rz(-2.984798) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.110072) q[1];
sx q[1];
rz(-1.8247461) q[1];
sx q[1];
rz(0.54540789) q[1];
rz(-pi) q[2];
rz(-2.429871) q[3];
sx q[3];
rz(-0.91152836) q[3];
sx q[3];
rz(-3.1275415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.86924187) q[2];
sx q[2];
rz(-1.7607471) q[2];
sx q[2];
rz(2.8273919) q[2];
rz(-2.3172486) q[3];
sx q[3];
rz(-0.45212513) q[3];
sx q[3];
rz(2.3468988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070351275) q[0];
sx q[0];
rz(-0.059878778) q[0];
sx q[0];
rz(-1.2605793) q[0];
rz(-2.4977327) q[1];
sx q[1];
rz(-1.9088129) q[1];
sx q[1];
rz(-0.018928122) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8986329) q[0];
sx q[0];
rz(-1.3071994) q[0];
sx q[0];
rz(1.5893057) q[0];
rz(-0.55982121) q[2];
sx q[2];
rz(-0.18351843) q[2];
sx q[2];
rz(0.57932094) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.85266528) q[1];
sx q[1];
rz(-1.3864494) q[1];
sx q[1];
rz(1.5194555) q[1];
rz(-1.8065679) q[3];
sx q[3];
rz(-1.8528432) q[3];
sx q[3];
rz(2.3896133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.54946047) q[2];
sx q[2];
rz(-0.38828725) q[2];
sx q[2];
rz(-2.4712759) q[2];
rz(2.629225) q[3];
sx q[3];
rz(-1.7497601) q[3];
sx q[3];
rz(1.7211154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55384127) q[0];
sx q[0];
rz(-1.8974263) q[0];
sx q[0];
rz(0.49945369) q[0];
rz(-1.5746501) q[1];
sx q[1];
rz(-0.27856871) q[1];
sx q[1];
rz(-2.0589028) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64718819) q[0];
sx q[0];
rz(-1.5218966) q[0];
sx q[0];
rz(2.5542269) q[0];
x q[1];
rz(0.033109025) q[2];
sx q[2];
rz(-0.81956714) q[2];
sx q[2];
rz(0.59226743) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.618256) q[1];
sx q[1];
rz(-0.54392951) q[1];
sx q[1];
rz(-3.1052599) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4389078) q[3];
sx q[3];
rz(-2.6583238) q[3];
sx q[3];
rz(1.0815222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7426976) q[2];
sx q[2];
rz(-1.9766786) q[2];
sx q[2];
rz(-2.6297074) q[2];
rz(2.7486457) q[3];
sx q[3];
rz(-1.7397375) q[3];
sx q[3];
rz(1.0569364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.186541) q[0];
sx q[0];
rz(-1.3165836) q[0];
sx q[0];
rz(-2.5008428) q[0];
rz(2.4004249) q[1];
sx q[1];
rz(-2.3186431) q[1];
sx q[1];
rz(2.9021312) q[1];
rz(-2.1869833) q[2];
sx q[2];
rz(-1.2381427) q[2];
sx q[2];
rz(2.8340813) q[2];
rz(-1.742733) q[3];
sx q[3];
rz(-0.83172432) q[3];
sx q[3];
rz(-1.5749501) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
