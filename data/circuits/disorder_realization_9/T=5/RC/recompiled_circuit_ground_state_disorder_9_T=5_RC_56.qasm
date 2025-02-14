OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.61152148) q[0];
sx q[0];
rz(-1.0728711) q[0];
sx q[0];
rz(1.7142417) q[0];
rz(-2.144835) q[1];
sx q[1];
rz(-0.61350322) q[1];
sx q[1];
rz(-0.065453425) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4213227) q[0];
sx q[0];
rz(-1.7593059) q[0];
sx q[0];
rz(-2.6557198) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5292306) q[2];
sx q[2];
rz(-1.9392559) q[2];
sx q[2];
rz(-0.17247904) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0411543) q[1];
sx q[1];
rz(-1.8722495) q[1];
sx q[1];
rz(0.90561066) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0569589) q[3];
sx q[3];
rz(-1.0118985) q[3];
sx q[3];
rz(-1.8098179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2685711) q[2];
sx q[2];
rz(-0.34653386) q[2];
sx q[2];
rz(-1.2968501) q[2];
rz(1.1067363) q[3];
sx q[3];
rz(-2.1745671) q[3];
sx q[3];
rz(1.6905674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1065555) q[0];
sx q[0];
rz(-2.4132001) q[0];
sx q[0];
rz(1.6075217) q[0];
rz(-0.84397498) q[1];
sx q[1];
rz(-1.5183828) q[1];
sx q[1];
rz(0.27110505) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26626884) q[0];
sx q[0];
rz(-1.3014587) q[0];
sx q[0];
rz(-1.6104417) q[0];
x q[1];
rz(1.9366802) q[2];
sx q[2];
rz(-0.26981631) q[2];
sx q[2];
rz(-1.3617977) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7162937) q[1];
sx q[1];
rz(-1.6669841) q[1];
sx q[1];
rz(2.7240818) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4976147) q[3];
sx q[3];
rz(-1.4212307) q[3];
sx q[3];
rz(0.73313802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.33807445) q[2];
sx q[2];
rz(-0.98793554) q[2];
sx q[2];
rz(0.79338497) q[2];
rz(0.042898305) q[3];
sx q[3];
rz(-1.2267313) q[3];
sx q[3];
rz(0.66810098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0778097) q[0];
sx q[0];
rz(-0.48352799) q[0];
sx q[0];
rz(-3.0384645) q[0];
rz(-2.8887796) q[1];
sx q[1];
rz(-1.4272855) q[1];
sx q[1];
rz(2.8275183) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8738587) q[0];
sx q[0];
rz(-1.9600057) q[0];
sx q[0];
rz(0.53797526) q[0];
rz(-1.6797596) q[2];
sx q[2];
rz(-0.48479947) q[2];
sx q[2];
rz(1.4919748) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9772226) q[1];
sx q[1];
rz(-0.80406351) q[1];
sx q[1];
rz(2.3146446) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57468412) q[3];
sx q[3];
rz(-2.8369224) q[3];
sx q[3];
rz(1.6449354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.126943) q[2];
sx q[2];
rz(-0.64121556) q[2];
sx q[2];
rz(-0.81652299) q[2];
rz(2.4294295) q[3];
sx q[3];
rz(-1.687259) q[3];
sx q[3];
rz(0.85675353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.6461058) q[0];
sx q[0];
rz(-3.0480338) q[0];
sx q[0];
rz(-2.5529472) q[0];
rz(0.37216392) q[1];
sx q[1];
rz(-1.2969505) q[1];
sx q[1];
rz(-2.1968496) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3550081) q[0];
sx q[0];
rz(-1.8898801) q[0];
sx q[0];
rz(0.068885013) q[0];
rz(2.0633374) q[2];
sx q[2];
rz(-1.8027935) q[2];
sx q[2];
rz(-1.7556695) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0081629) q[1];
sx q[1];
rz(-1.1436698) q[1];
sx q[1];
rz(0.62469805) q[1];
rz(-pi) q[2];
rz(-2.9507941) q[3];
sx q[3];
rz(-1.2114297) q[3];
sx q[3];
rz(2.1426853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9113691) q[2];
sx q[2];
rz(-1.6065803) q[2];
sx q[2];
rz(1.6564507) q[2];
rz(2.7459775) q[3];
sx q[3];
rz(-1.6499358) q[3];
sx q[3];
rz(0.0609456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38451251) q[0];
sx q[0];
rz(-2.7482996) q[0];
sx q[0];
rz(-1.242189) q[0];
rz(0.077266129) q[1];
sx q[1];
rz(-1.2178414) q[1];
sx q[1];
rz(3.0526615) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0384232) q[0];
sx q[0];
rz(-2.9154202) q[0];
sx q[0];
rz(0.75407501) q[0];
rz(-pi) q[1];
rz(-1.4706663) q[2];
sx q[2];
rz(-1.0271974) q[2];
sx q[2];
rz(1.3058654) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4684852) q[1];
sx q[1];
rz(-1.2947646) q[1];
sx q[1];
rz(-2.0620146) q[1];
rz(3.0432391) q[3];
sx q[3];
rz(-1.438201) q[3];
sx q[3];
rz(-1.7817532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5605805) q[2];
sx q[2];
rz(-1.7580914) q[2];
sx q[2];
rz(1.2233454) q[2];
rz(1.9129725) q[3];
sx q[3];
rz(-2.2556428) q[3];
sx q[3];
rz(-0.74738735) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0864047) q[0];
sx q[0];
rz(-2.3265657) q[0];
sx q[0];
rz(-2.1888457) q[0];
rz(-1.3505666) q[1];
sx q[1];
rz(-1.9074651) q[1];
sx q[1];
rz(-1.9780673) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19630963) q[0];
sx q[0];
rz(-0.74271129) q[0];
sx q[0];
rz(-0.16603827) q[0];
rz(-pi) q[1];
rz(1.8356531) q[2];
sx q[2];
rz(-0.9914757) q[2];
sx q[2];
rz(0.91897041) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4965373) q[1];
sx q[1];
rz(-2.8926507) q[1];
sx q[1];
rz(1.7824936) q[1];
rz(2.5361193) q[3];
sx q[3];
rz(-1.2126727) q[3];
sx q[3];
rz(-0.78152657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20840883) q[2];
sx q[2];
rz(-1.9219425) q[2];
sx q[2];
rz(-1.5737994) q[2];
rz(2.9257704) q[3];
sx q[3];
rz(-0.62041557) q[3];
sx q[3];
rz(-0.1568493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
rz(-0.02483524) q[0];
sx q[0];
rz(-0.76501608) q[0];
sx q[0];
rz(-1.5417954) q[0];
rz(2.722591) q[1];
sx q[1];
rz(-1.1083138) q[1];
sx q[1];
rz(1.3655183) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1772108) q[0];
sx q[0];
rz(-2.1809077) q[0];
sx q[0];
rz(-2.9047853) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0232685) q[2];
sx q[2];
rz(-0.86811262) q[2];
sx q[2];
rz(-1.4066525) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1394386) q[1];
sx q[1];
rz(-1.4055058) q[1];
sx q[1];
rz(-2.2846245) q[1];
x q[2];
rz(0.98181866) q[3];
sx q[3];
rz(-1.8659628) q[3];
sx q[3];
rz(1.3491614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8968203) q[2];
sx q[2];
rz(-0.12910566) q[2];
sx q[2];
rz(-1.6884035) q[2];
rz(1.3990654) q[3];
sx q[3];
rz(-1.1996256) q[3];
sx q[3];
rz(2.3340268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6145265) q[0];
sx q[0];
rz(-0.87166059) q[0];
sx q[0];
rz(0.51323071) q[0];
rz(0.97663438) q[1];
sx q[1];
rz(-2.4726424) q[1];
sx q[1];
rz(1.4248779) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34951613) q[0];
sx q[0];
rz(-2.2648025) q[0];
sx q[0];
rz(2.9727139) q[0];
rz(-pi) q[1];
rz(0.11995671) q[2];
sx q[2];
rz(-0.55177125) q[2];
sx q[2];
rz(0.063869501) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1816563) q[1];
sx q[1];
rz(-1.8605369) q[1];
sx q[1];
rz(2.0253889) q[1];
rz(-3.023671) q[3];
sx q[3];
rz(-1.4050845) q[3];
sx q[3];
rz(2.5664701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.49154115) q[2];
sx q[2];
rz(-2.1341133) q[2];
sx q[2];
rz(0.36637351) q[2];
rz(-2.5036687) q[3];
sx q[3];
rz(-1.7584691) q[3];
sx q[3];
rz(-1.6320419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3731308) q[0];
sx q[0];
rz(-0.32485425) q[0];
sx q[0];
rz(-1.984206) q[0];
rz(-2.6034082) q[1];
sx q[1];
rz(-1.3342074) q[1];
sx q[1];
rz(-0.023712637) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0932281) q[0];
sx q[0];
rz(-0.95010883) q[0];
sx q[0];
rz(-1.4976504) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41466434) q[2];
sx q[2];
rz(-0.86277308) q[2];
sx q[2];
rz(1.9239349) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7754375) q[1];
sx q[1];
rz(-1.0682271) q[1];
sx q[1];
rz(1.9910046) q[1];
x q[2];
rz(-3.0160041) q[3];
sx q[3];
rz(-1.3712582) q[3];
sx q[3];
rz(-0.34717595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1553161) q[2];
sx q[2];
rz(-1.3331058) q[2];
sx q[2];
rz(-0.39815608) q[2];
rz(2.6768173) q[3];
sx q[3];
rz(-2.0538752) q[3];
sx q[3];
rz(1.8006511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-1.583113) q[0];
sx q[0];
rz(-0.3322424) q[0];
sx q[0];
rz(2.1930021) q[0];
rz(0.96725431) q[1];
sx q[1];
rz(-1.1525258) q[1];
sx q[1];
rz(2.1327877) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7083884) q[0];
sx q[0];
rz(-2.3713667) q[0];
sx q[0];
rz(1.5807093) q[0];
rz(-2.45935) q[2];
sx q[2];
rz(-2.3385404) q[2];
sx q[2];
rz(2.19953) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0151143) q[1];
sx q[1];
rz(-0.51869828) q[1];
sx q[1];
rz(0.37094231) q[1];
rz(-pi) q[2];
x q[2];
rz(2.860805) q[3];
sx q[3];
rz(-0.6774582) q[3];
sx q[3];
rz(1.2408797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41445109) q[2];
sx q[2];
rz(-2.264617) q[2];
sx q[2];
rz(0.021075185) q[2];
rz(2.6575139) q[3];
sx q[3];
rz(-2.0120072) q[3];
sx q[3];
rz(0.67897236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2912343) q[0];
sx q[0];
rz(-1.5342916) q[0];
sx q[0];
rz(1.7271484) q[0];
rz(0.5961295) q[1];
sx q[1];
rz(-2.3565751) q[1];
sx q[1];
rz(-0.48859488) q[1];
rz(2.1837744) q[2];
sx q[2];
rz(-0.25479813) q[2];
sx q[2];
rz(0.23164498) q[2];
rz(-3.0335043) q[3];
sx q[3];
rz(-1.3664403) q[3];
sx q[3];
rz(-1.4667778) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
