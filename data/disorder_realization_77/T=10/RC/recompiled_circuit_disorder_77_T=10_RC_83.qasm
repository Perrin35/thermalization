OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(-0.84364426) q[0];
sx q[0];
rz(0.16790976) q[0];
rz(-1.9703938) q[1];
sx q[1];
rz(-0.29532239) q[1];
sx q[1];
rz(3.0854316) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0060318) q[0];
sx q[0];
rz(-0.52838415) q[0];
sx q[0];
rz(-2.0023268) q[0];
rz(-0.21284717) q[2];
sx q[2];
rz(-2.2058862) q[2];
sx q[2];
rz(2.0095306) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.010477) q[1];
sx q[1];
rz(-1.3033723) q[1];
sx q[1];
rz(2.1285776) q[1];
rz(-0.16529103) q[3];
sx q[3];
rz(-2.878302) q[3];
sx q[3];
rz(-0.96112448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7636259) q[2];
sx q[2];
rz(-0.28181919) q[2];
sx q[2];
rz(0.4326694) q[2];
rz(1.9487322) q[3];
sx q[3];
rz(-1.9038707) q[3];
sx q[3];
rz(-2.7584934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52779657) q[0];
sx q[0];
rz(-0.48848099) q[0];
sx q[0];
rz(-1.8288076) q[0];
rz(-2.9361172) q[1];
sx q[1];
rz(-0.97646362) q[1];
sx q[1];
rz(1.9899433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8173556) q[0];
sx q[0];
rz(-1.865987) q[0];
sx q[0];
rz(-2.2833707) q[0];
rz(-pi) q[1];
rz(0.58354124) q[2];
sx q[2];
rz(-1.1569287) q[2];
sx q[2];
rz(1.5577424) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.68576751) q[1];
sx q[1];
rz(-1.1357422) q[1];
sx q[1];
rz(1.1063834) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50206708) q[3];
sx q[3];
rz(-1.2289398) q[3];
sx q[3];
rz(1.348192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1318704) q[2];
sx q[2];
rz(-1.4902318) q[2];
sx q[2];
rz(-2.9197664) q[2];
rz(0.37718537) q[3];
sx q[3];
rz(-0.42755869) q[3];
sx q[3];
rz(0.77243531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8310228) q[0];
sx q[0];
rz(-0.092386827) q[0];
sx q[0];
rz(-0.036852766) q[0];
rz(-0.82551461) q[1];
sx q[1];
rz(-1.8258391) q[1];
sx q[1];
rz(-0.056578606) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.409257) q[0];
sx q[0];
rz(-1.1261254) q[0];
sx q[0];
rz(2.1952941) q[0];
rz(-0.47358863) q[2];
sx q[2];
rz(-2.8550672) q[2];
sx q[2];
rz(-2.5055656) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0128855) q[1];
sx q[1];
rz(-1.371908) q[1];
sx q[1];
rz(0.30602869) q[1];
rz(-pi) q[2];
rz(0.49719663) q[3];
sx q[3];
rz(-0.70780863) q[3];
sx q[3];
rz(-1.4149815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.27292192) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(0.92612129) q[2];
rz(-0.55666322) q[3];
sx q[3];
rz(-0.29354468) q[3];
sx q[3];
rz(-2.0986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91519231) q[0];
sx q[0];
rz(-2.4214348) q[0];
sx q[0];
rz(-2.3994989) q[0];
rz(2.0023951) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(0.46359584) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0915506) q[0];
sx q[0];
rz(-0.76489641) q[0];
sx q[0];
rz(1.1372304) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8088403) q[2];
sx q[2];
rz(-1.5126265) q[2];
sx q[2];
rz(1.0156877) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1222893) q[1];
sx q[1];
rz(-2.3483854) q[1];
sx q[1];
rz(1.8122458) q[1];
rz(2.6632518) q[3];
sx q[3];
rz(-1.0676749) q[3];
sx q[3];
rz(-0.92418811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3670369) q[2];
sx q[2];
rz(-2.98428) q[2];
sx q[2];
rz(3.0920933) q[2];
rz(3.0130623) q[3];
sx q[3];
rz(-1.5481719) q[3];
sx q[3];
rz(-3.0226504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0054935) q[0];
sx q[0];
rz(-2.7042784) q[0];
sx q[0];
rz(0.29770011) q[0];
rz(-2.659335) q[1];
sx q[1];
rz(-2.3869956) q[1];
sx q[1];
rz(-0.94435) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0194861) q[0];
sx q[0];
rz(-1.5257611) q[0];
sx q[0];
rz(-1.7800063) q[0];
rz(-pi) q[1];
x q[1];
rz(1.96825) q[2];
sx q[2];
rz(-2.3059418) q[2];
sx q[2];
rz(0.022692516) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.797232) q[1];
sx q[1];
rz(-1.6318983) q[1];
sx q[1];
rz(-3.1394715) q[1];
rz(-2.7086908) q[3];
sx q[3];
rz(-1.8042943) q[3];
sx q[3];
rz(1.7683065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8828316) q[2];
sx q[2];
rz(-2.0488887) q[2];
sx q[2];
rz(-3.0333701) q[2];
rz(0.0023068874) q[3];
sx q[3];
rz(-1.5294411) q[3];
sx q[3];
rz(2.8172857) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7047983) q[0];
sx q[0];
rz(-0.44610766) q[0];
sx q[0];
rz(-0.58445245) q[0];
rz(-2.2553518) q[1];
sx q[1];
rz(-2.5247572) q[1];
sx q[1];
rz(-3.086673) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14558218) q[0];
sx q[0];
rz(-1.8390363) q[0];
sx q[0];
rz(-0.085573816) q[0];
x q[1];
rz(0.018718406) q[2];
sx q[2];
rz(-1.2476377) q[2];
sx q[2];
rz(-1.2013555) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1153032) q[1];
sx q[1];
rz(-0.68740986) q[1];
sx q[1];
rz(0.61647146) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46338007) q[3];
sx q[3];
rz(-1.0372835) q[3];
sx q[3];
rz(-2.2802071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3871258) q[2];
sx q[2];
rz(-3.0209164) q[2];
sx q[2];
rz(-2.1248655) q[2];
rz(-0.54404849) q[3];
sx q[3];
rz(-0.35990158) q[3];
sx q[3];
rz(-1.8090766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6761557) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(-2.8570535) q[0];
rz(-0.94447213) q[1];
sx q[1];
rz(-1.9453134) q[1];
sx q[1];
rz(2.231266) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9757662) q[0];
sx q[0];
rz(-1.5162139) q[0];
sx q[0];
rz(1.5058917) q[0];
rz(-pi) q[1];
rz(-0.89332135) q[2];
sx q[2];
rz(-2.6258694) q[2];
sx q[2];
rz(-0.90781462) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6701723) q[1];
sx q[1];
rz(-0.52392611) q[1];
sx q[1];
rz(0.36797221) q[1];
rz(-pi) q[2];
rz(-2.4207553) q[3];
sx q[3];
rz(-1.1158873) q[3];
sx q[3];
rz(0.37973675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7515144) q[2];
sx q[2];
rz(-0.063515924) q[2];
sx q[2];
rz(-0.92203036) q[2];
rz(-2.5743124) q[3];
sx q[3];
rz(-1.7276948) q[3];
sx q[3];
rz(-2.1218307) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89408016) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(3.0122053) q[0];
rz(2.5091876) q[1];
sx q[1];
rz(-2.1148732) q[1];
sx q[1];
rz(2.8410889) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6209517) q[0];
sx q[0];
rz(-1.9410987) q[0];
sx q[0];
rz(-2.305549) q[0];
rz(-pi) q[1];
rz(0.76086107) q[2];
sx q[2];
rz(-1.0957452) q[2];
sx q[2];
rz(2.8617815) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.955519) q[1];
sx q[1];
rz(-1.7454073) q[1];
sx q[1];
rz(2.7956635) q[1];
x q[2];
rz(-1.9037876) q[3];
sx q[3];
rz(-0.78562842) q[3];
sx q[3];
rz(-1.1803407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5552716) q[2];
sx q[2];
rz(-0.99493146) q[2];
sx q[2];
rz(-0.78197455) q[2];
rz(-0.55109763) q[3];
sx q[3];
rz(-1.3827773) q[3];
sx q[3];
rz(2.9836695) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5732116) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(3.0138299) q[0];
rz(2.5993775) q[1];
sx q[1];
rz(-2.1844889) q[1];
sx q[1];
rz(2.382747) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.471506) q[0];
sx q[0];
rz(-0.42450464) q[0];
sx q[0];
rz(1.5295117) q[0];
x q[1];
rz(1.3779638) q[2];
sx q[2];
rz(-0.97059965) q[2];
sx q[2];
rz(2.6892975) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.89275937) q[1];
sx q[1];
rz(-1.7264688) q[1];
sx q[1];
rz(0.70181904) q[1];
rz(0.014563668) q[3];
sx q[3];
rz(-2.2363538) q[3];
sx q[3];
rz(-1.2138106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1252497) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(0.49003595) q[2];
rz(1.7193433) q[3];
sx q[3];
rz(-1.2079206) q[3];
sx q[3];
rz(-2.0786044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7816417) q[0];
sx q[0];
rz(-0.61976969) q[0];
sx q[0];
rz(0.075335659) q[0];
rz(-0.8967337) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(-2.5316701) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2595554) q[0];
sx q[0];
rz(-1.2872739) q[0];
sx q[0];
rz(-2.3638704) q[0];
x q[1];
rz(-1.8462734) q[2];
sx q[2];
rz(-1.2071949) q[2];
sx q[2];
rz(-1.7737349) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3366821) q[1];
sx q[1];
rz(-1.8098127) q[1];
sx q[1];
rz(-2.3906624) q[1];
rz(-pi) q[2];
rz(-2.6033953) q[3];
sx q[3];
rz(-0.81848577) q[3];
sx q[3];
rz(-1.8173816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.909409) q[2];
sx q[2];
rz(-0.82513088) q[2];
sx q[2];
rz(-0.71371901) q[2];
rz(0.37832007) q[3];
sx q[3];
rz(-2.648073) q[3];
sx q[3];
rz(-2.2617214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80355766) q[0];
sx q[0];
rz(-1.150158) q[0];
sx q[0];
rz(-1.5858571) q[0];
rz(-0.65080416) q[1];
sx q[1];
rz(-1.4918068) q[1];
sx q[1];
rz(3.0202958) q[1];
rz(-3.1109839) q[2];
sx q[2];
rz(-1.7666713) q[2];
sx q[2];
rz(-0.91790744) q[2];
rz(1.0999023) q[3];
sx q[3];
rz(-1.3445911) q[3];
sx q[3];
rz(0.54683987) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
