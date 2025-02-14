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
rz(2.0119014) q[0];
rz(-2.5870485) q[1];
sx q[1];
rz(-2.6362042) q[1];
sx q[1];
rz(0.55355054) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3239297) q[0];
sx q[0];
rz(-0.50369064) q[0];
sx q[0];
rz(-1.5970741) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0099382691) q[2];
sx q[2];
rz(-1.9807837) q[2];
sx q[2];
rz(2.8824497) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5207102) q[1];
sx q[1];
rz(-1.1739879) q[1];
sx q[1];
rz(2.3842393) q[1];
rz(-pi) q[2];
rz(-0.042993892) q[3];
sx q[3];
rz(-1.3519796) q[3];
sx q[3];
rz(2.9584663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1393789) q[2];
sx q[2];
rz(-1.4993818) q[2];
sx q[2];
rz(1.9726944) q[2];
rz(-0.72541952) q[3];
sx q[3];
rz(-1.231671) q[3];
sx q[3];
rz(1.5861082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.91189522) q[0];
sx q[0];
rz(-2.732369) q[0];
sx q[0];
rz(-1.0503861) q[0];
rz(-1.3788252) q[1];
sx q[1];
rz(-1.4070815) q[1];
sx q[1];
rz(1.001766) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40658298) q[0];
sx q[0];
rz(-1.4676784) q[0];
sx q[0];
rz(1.3117657) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0323192) q[2];
sx q[2];
rz(-1.8428728) q[2];
sx q[2];
rz(2.9124454) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.54047062) q[1];
sx q[1];
rz(-1.8007638) q[1];
sx q[1];
rz(2.7011987) q[1];
rz(-pi) q[2];
rz(2.7185387) q[3];
sx q[3];
rz(-1.3025369) q[3];
sx q[3];
rz(-0.34706193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0687781) q[2];
sx q[2];
rz(-0.38663703) q[2];
sx q[2];
rz(-1.4978503) q[2];
rz(2.1554598) q[3];
sx q[3];
rz(-2.17406) q[3];
sx q[3];
rz(0.51924527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3471974) q[0];
sx q[0];
rz(-0.47413844) q[0];
sx q[0];
rz(-0.073632181) q[0];
rz(2.4108048) q[1];
sx q[1];
rz(-1.7981139) q[1];
sx q[1];
rz(0.98339072) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7950566) q[0];
sx q[0];
rz(-1.4162401) q[0];
sx q[0];
rz(-1.2861504) q[0];
x q[1];
rz(2.2762301) q[2];
sx q[2];
rz(-1.1956805) q[2];
sx q[2];
rz(-2.4789711) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4913339) q[1];
sx q[1];
rz(-2.8425807) q[1];
sx q[1];
rz(2.0102802) q[1];
rz(-pi) q[2];
rz(2.0737209) q[3];
sx q[3];
rz(-1.6391603) q[3];
sx q[3];
rz(1.158744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.44707766) q[2];
sx q[2];
rz(-2.647001) q[2];
sx q[2];
rz(-0.05973235) q[2];
rz(-2.6308502) q[3];
sx q[3];
rz(-2.7849438) q[3];
sx q[3];
rz(1.7662778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.315595) q[0];
sx q[0];
rz(-2.4331278) q[0];
sx q[0];
rz(-2.3009756) q[0];
rz(-2.9261342) q[1];
sx q[1];
rz(-2.0550199) q[1];
sx q[1];
rz(-3.0256909) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8944696) q[0];
sx q[0];
rz(-0.43465675) q[0];
sx q[0];
rz(0.91095319) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5253549) q[2];
sx q[2];
rz(-0.35348693) q[2];
sx q[2];
rz(2.7797159) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2167989) q[1];
sx q[1];
rz(-1.8571641) q[1];
sx q[1];
rz(0.91735265) q[1];
x q[2];
rz(-0.23502879) q[3];
sx q[3];
rz(-2.7209366) q[3];
sx q[3];
rz(3.0893081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0243715) q[2];
sx q[2];
rz(-1.8481959) q[2];
sx q[2];
rz(0.8563861) q[2];
rz(-0.74445009) q[3];
sx q[3];
rz(-2.2529633) q[3];
sx q[3];
rz(0.3652679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(0.15904467) q[0];
sx q[0];
rz(-2.5578975) q[0];
sx q[0];
rz(2.1481376) q[0];
rz(-1.1194725) q[1];
sx q[1];
rz(-1.9529724) q[1];
sx q[1];
rz(-3.1192034) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3684621) q[0];
sx q[0];
rz(-1.4903632) q[0];
sx q[0];
rz(0.12534938) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3359433) q[2];
sx q[2];
rz(-2.4107526) q[2];
sx q[2];
rz(-1.8455055) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4926947) q[1];
sx q[1];
rz(-1.9580868) q[1];
sx q[1];
rz(2.5084916) q[1];
rz(-0.56110142) q[3];
sx q[3];
rz(-0.77563876) q[3];
sx q[3];
rz(2.5833094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8523954) q[2];
sx q[2];
rz(-2.9183233) q[2];
sx q[2];
rz(-1.4168463) q[2];
rz(1.019574) q[3];
sx q[3];
rz(-0.9905147) q[3];
sx q[3];
rz(0.31653252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1932909) q[0];
sx q[0];
rz(-0.010007771) q[0];
sx q[0];
rz(2.365812) q[0];
rz(-1.4099482) q[1];
sx q[1];
rz(-1.3103139) q[1];
sx q[1];
rz(1.5820679) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5821302) q[0];
sx q[0];
rz(-3.0728615) q[0];
sx q[0];
rz(1.1148081) q[0];
rz(-pi) q[1];
rz(-0.36119097) q[2];
sx q[2];
rz(-1.0345999) q[2];
sx q[2];
rz(2.9597832) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.74461339) q[1];
sx q[1];
rz(-1.1049264) q[1];
sx q[1];
rz(-0.14772285) q[1];
rz(-1.9385064) q[3];
sx q[3];
rz(-2.0597405) q[3];
sx q[3];
rz(2.6206102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55283028) q[2];
sx q[2];
rz(-1.7071807) q[2];
sx q[2];
rz(2.8130048) q[2];
rz(1.5821247) q[3];
sx q[3];
rz(-0.72454536) q[3];
sx q[3];
rz(-2.4174378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62411031) q[0];
sx q[0];
rz(-0.8172074) q[0];
sx q[0];
rz(0.66993129) q[0];
rz(-1.5176557) q[1];
sx q[1];
rz(-2.0567963) q[1];
sx q[1];
rz(-2.0673015) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5094604) q[0];
sx q[0];
rz(-1.7863417) q[0];
sx q[0];
rz(2.4556405) q[0];
rz(-pi) q[1];
rz(-2.7781717) q[2];
sx q[2];
rz(-1.7426024) q[2];
sx q[2];
rz(1.8653143) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8855369) q[1];
sx q[1];
rz(-0.32840751) q[1];
sx q[1];
rz(-1.8610015) q[1];
x q[2];
rz(-1.8378476) q[3];
sx q[3];
rz(-0.69483583) q[3];
sx q[3];
rz(1.9766493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7488148) q[2];
sx q[2];
rz(-0.98069507) q[2];
sx q[2];
rz(1.7559715) q[2];
rz(-2.4393926) q[3];
sx q[3];
rz(-0.32679138) q[3];
sx q[3];
rz(2.0328111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4362519) q[0];
sx q[0];
rz(-1.5974644) q[0];
sx q[0];
rz(-2.3633549) q[0];
rz(1.8477919) q[1];
sx q[1];
rz(-1.2575282) q[1];
sx q[1];
rz(3.0857118) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27625448) q[0];
sx q[0];
rz(-2.8408065) q[0];
sx q[0];
rz(0.24888034) q[0];
rz(-pi) q[1];
rz(-2.915809) q[2];
sx q[2];
rz(-1.3874545) q[2];
sx q[2];
rz(-2.8933337) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3371007) q[1];
sx q[1];
rz(-2.7332515) q[1];
sx q[1];
rz(0.33894811) q[1];
x q[2];
rz(2.6229739) q[3];
sx q[3];
rz(-2.2317776) q[3];
sx q[3];
rz(2.3346614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0539661) q[2];
sx q[2];
rz(-2.1199333) q[2];
sx q[2];
rz(2.450024) q[2];
rz(-2.9223053) q[3];
sx q[3];
rz(-1.4838452) q[3];
sx q[3];
rz(-1.9482435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9241065) q[0];
sx q[0];
rz(-0.036343887) q[0];
sx q[0];
rz(2.0935667) q[0];
rz(0.74608392) q[1];
sx q[1];
rz(-1.8522976) q[1];
sx q[1];
rz(0.41864905) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11514535) q[0];
sx q[0];
rz(-1.5799684) q[0];
sx q[0];
rz(1.5546726) q[0];
rz(-0.23048909) q[2];
sx q[2];
rz(-0.93991342) q[2];
sx q[2];
rz(2.0216551) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.56142226) q[1];
sx q[1];
rz(-0.19374312) q[1];
sx q[1];
rz(-2.3592739) q[1];
rz(-pi) q[2];
rz(-3.0935086) q[3];
sx q[3];
rz(-0.86838956) q[3];
sx q[3];
rz(0.66437214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5082385) q[2];
sx q[2];
rz(-0.60428047) q[2];
sx q[2];
rz(0.61887997) q[2];
rz(-2.6887756) q[3];
sx q[3];
rz(-1.489095) q[3];
sx q[3];
rz(0.49310163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(2.598269) q[0];
sx q[0];
rz(-3.0040574) q[0];
sx q[0];
rz(1.2602873) q[0];
rz(-2.8434985) q[1];
sx q[1];
rz(-0.9577848) q[1];
sx q[1];
rz(1.6424929) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84628203) q[0];
sx q[0];
rz(-0.75937004) q[0];
sx q[0];
rz(-1.6037386) q[0];
x q[1];
rz(0.46447931) q[2];
sx q[2];
rz(-0.3128007) q[2];
sx q[2];
rz(-2.0020747) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.11718845) q[1];
sx q[1];
rz(-2.9147101) q[1];
sx q[1];
rz(-1.6443166) q[1];
x q[2];
rz(-1.1785678) q[3];
sx q[3];
rz(-1.7081327) q[3];
sx q[3];
rz(1.3016421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5160211) q[2];
sx q[2];
rz(-1.8578153) q[2];
sx q[2];
rz(0.18216356) q[2];
rz(-0.22570172) q[3];
sx q[3];
rz(-1.00939) q[3];
sx q[3];
rz(-1.3753447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7094649) q[0];
sx q[0];
rz(-1.3668677) q[0];
sx q[0];
rz(0.9912542) q[0];
rz(-1.9215012) q[1];
sx q[1];
rz(-2.5972001) q[1];
sx q[1];
rz(-0.77144815) q[1];
rz(1.3924072) q[2];
sx q[2];
rz(-1.6606944) q[2];
sx q[2];
rz(0.14002945) q[2];
rz(-0.75481244) q[3];
sx q[3];
rz(-0.94554286) q[3];
sx q[3];
rz(1.7331358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
