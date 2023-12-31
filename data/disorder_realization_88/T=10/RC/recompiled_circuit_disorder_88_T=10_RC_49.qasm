OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8945583) q[0];
sx q[0];
rz(-1.2556827) q[0];
sx q[0];
rz(-0.32796252) q[0];
rz(2.9070931) q[1];
sx q[1];
rz(-0.20107888) q[1];
sx q[1];
rz(0.091436401) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8592005) q[0];
sx q[0];
rz(-0.99635591) q[0];
sx q[0];
rz(-2.0496856) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.882471) q[2];
sx q[2];
rz(-1.4403575) q[2];
sx q[2];
rz(0.63333095) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2135703) q[1];
sx q[1];
rz(-1.808128) q[1];
sx q[1];
rz(3.0568394) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9760194) q[3];
sx q[3];
rz(-1.3888437) q[3];
sx q[3];
rz(2.6973157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4771007) q[2];
sx q[2];
rz(-0.97350073) q[2];
sx q[2];
rz(-1.1260024) q[2];
rz(0.27515718) q[3];
sx q[3];
rz(-2.5313009) q[3];
sx q[3];
rz(2.2385105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0098410957) q[0];
sx q[0];
rz(-0.69350243) q[0];
sx q[0];
rz(0.69357187) q[0];
rz(1.0961078) q[1];
sx q[1];
rz(-0.98384905) q[1];
sx q[1];
rz(0.19031659) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099146518) q[0];
sx q[0];
rz(-2.6876246) q[0];
sx q[0];
rz(-2.0367665) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6066577) q[2];
sx q[2];
rz(-1.2282279) q[2];
sx q[2];
rz(-1.5869706) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39365921) q[1];
sx q[1];
rz(-0.84888443) q[1];
sx q[1];
rz(0.02971239) q[1];
rz(-2.4137647) q[3];
sx q[3];
rz(-1.9350633) q[3];
sx q[3];
rz(1.4621853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6341614) q[2];
sx q[2];
rz(-0.71175152) q[2];
sx q[2];
rz(1.2197536) q[2];
rz(-2.9988585) q[3];
sx q[3];
rz(-1.0120564) q[3];
sx q[3];
rz(-2.232961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74293566) q[0];
sx q[0];
rz(-1.0899028) q[0];
sx q[0];
rz(0.8202585) q[0];
rz(-0.29207686) q[1];
sx q[1];
rz(-2.0675817) q[1];
sx q[1];
rz(1.2480199) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4688063) q[0];
sx q[0];
rz(-0.60657036) q[0];
sx q[0];
rz(0.16288217) q[0];
x q[1];
rz(1.5065932) q[2];
sx q[2];
rz(-1.3639796) q[2];
sx q[2];
rz(-1.9931672) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.65590811) q[1];
sx q[1];
rz(-0.67042065) q[1];
sx q[1];
rz(0.1078492) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2218024) q[3];
sx q[3];
rz(-1.8120822) q[3];
sx q[3];
rz(1.9581219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32039207) q[2];
sx q[2];
rz(-2.0911066) q[2];
sx q[2];
rz(-1.9281663) q[2];
rz(-2.976867) q[3];
sx q[3];
rz(-0.89403331) q[3];
sx q[3];
rz(-0.081610672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7320025) q[0];
sx q[0];
rz(-1.859917) q[0];
sx q[0];
rz(-0.18606342) q[0];
rz(2.9371254) q[1];
sx q[1];
rz(-2.6696413) q[1];
sx q[1];
rz(-1.2971372) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38940934) q[0];
sx q[0];
rz(-1.7011189) q[0];
sx q[0];
rz(-0.043686314) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16764955) q[2];
sx q[2];
rz(-1.8561346) q[2];
sx q[2];
rz(-0.97380762) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.12249494) q[1];
sx q[1];
rz(-0.88224733) q[1];
sx q[1];
rz(2.3282939) q[1];
rz(-pi) q[2];
rz(-1.4947055) q[3];
sx q[3];
rz(-0.18008672) q[3];
sx q[3];
rz(2.0538581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2356448) q[2];
sx q[2];
rz(-0.79259688) q[2];
sx q[2];
rz(2.2223991) q[2];
rz(-0.32133189) q[3];
sx q[3];
rz(-1.0682169) q[3];
sx q[3];
rz(-1.8937768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.17386757) q[0];
sx q[0];
rz(-0.50650948) q[0];
sx q[0];
rz(0.87819535) q[0];
rz(-1.8163266) q[1];
sx q[1];
rz(-1.4122496) q[1];
sx q[1];
rz(-1.7153046) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.722055) q[0];
sx q[0];
rz(-0.66156045) q[0];
sx q[0];
rz(1.4369591) q[0];
rz(1.2582785) q[2];
sx q[2];
rz(-1.7973571) q[2];
sx q[2];
rz(0.3414008) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.95477415) q[1];
sx q[1];
rz(-2.2983119) q[1];
sx q[1];
rz(-2.9963042) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32894965) q[3];
sx q[3];
rz(-0.79862404) q[3];
sx q[3];
rz(-0.22508276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2720126) q[2];
sx q[2];
rz(-11*pi/13) q[2];
sx q[2];
rz(-0.041794725) q[2];
rz(-3.0801008) q[3];
sx q[3];
rz(-1.0890591) q[3];
sx q[3];
rz(-0.4367691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3549266) q[0];
sx q[0];
rz(-0.085952856) q[0];
sx q[0];
rz(2.1110995) q[0];
rz(-2.4018535) q[1];
sx q[1];
rz(-1.5286427) q[1];
sx q[1];
rz(-0.57156634) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7874998) q[0];
sx q[0];
rz(-2.0797074) q[0];
sx q[0];
rz(2.4718168) q[0];
rz(-1.7320485) q[2];
sx q[2];
rz(-0.67220062) q[2];
sx q[2];
rz(-1.6598998) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.9532277) q[1];
sx q[1];
rz(-1.1168915) q[1];
sx q[1];
rz(-2.8737349) q[1];
rz(-pi) q[2];
rz(-2.4404293) q[3];
sx q[3];
rz(-0.94651604) q[3];
sx q[3];
rz(-3.0091156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.968154) q[2];
sx q[2];
rz(-2.3996694) q[2];
sx q[2];
rz(2.5773933) q[2];
rz(-0.12600222) q[3];
sx q[3];
rz(-1.4583476) q[3];
sx q[3];
rz(-0.29461598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.0512222) q[0];
sx q[0];
rz(-2.9416961) q[0];
sx q[0];
rz(-2.4293161) q[0];
rz(0.5258711) q[1];
sx q[1];
rz(-0.41627517) q[1];
sx q[1];
rz(-2.4760822) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80488801) q[0];
sx q[0];
rz(-0.24222736) q[0];
sx q[0];
rz(1.5750984) q[0];
rz(-2.268157) q[2];
sx q[2];
rz(-1.4566112) q[2];
sx q[2];
rz(2.9261677) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4588786) q[1];
sx q[1];
rz(-2.6280118) q[1];
sx q[1];
rz(-2.0290124) q[1];
rz(1.8129559) q[3];
sx q[3];
rz(-1.7112964) q[3];
sx q[3];
rz(-1.5797918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9843288) q[2];
sx q[2];
rz(-0.66005808) q[2];
sx q[2];
rz(1.4228014) q[2];
rz(-0.11519365) q[3];
sx q[3];
rz(-2.6264103) q[3];
sx q[3];
rz(-2.9490024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049906235) q[0];
sx q[0];
rz(-2.005907) q[0];
sx q[0];
rz(2.9507622) q[0];
rz(-2.514839) q[1];
sx q[1];
rz(-1.0160867) q[1];
sx q[1];
rz(-0.33871067) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20913798) q[0];
sx q[0];
rz(-1.7212241) q[0];
sx q[0];
rz(1.9183137) q[0];
rz(-pi) q[1];
rz(0.4523925) q[2];
sx q[2];
rz(-2.6466742) q[2];
sx q[2];
rz(0.53168833) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.72570669) q[1];
sx q[1];
rz(-2.6542414) q[1];
sx q[1];
rz(1.332167) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8975774) q[3];
sx q[3];
rz(-1.3076926) q[3];
sx q[3];
rz(-2.9941878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.64615858) q[2];
sx q[2];
rz(-0.63445264) q[2];
sx q[2];
rz(-0.70043606) q[2];
rz(-2.2436079) q[3];
sx q[3];
rz(-1.286641) q[3];
sx q[3];
rz(2.9680796) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53196466) q[0];
sx q[0];
rz(-0.48159972) q[0];
sx q[0];
rz(2.4832446) q[0];
rz(0.61093962) q[1];
sx q[1];
rz(-1.2555723) q[1];
sx q[1];
rz(3.0019965) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4021554) q[0];
sx q[0];
rz(-1.6318775) q[0];
sx q[0];
rz(3.1215467) q[0];
rz(-pi) q[1];
rz(2.326968) q[2];
sx q[2];
rz(-1.5280208) q[2];
sx q[2];
rz(-0.2864366) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.86233625) q[1];
sx q[1];
rz(-1.5039872) q[1];
sx q[1];
rz(0.52211296) q[1];
x q[2];
rz(2.2152882) q[3];
sx q[3];
rz(-1.5288121) q[3];
sx q[3];
rz(-0.55461649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5247941) q[2];
sx q[2];
rz(-2.1308265) q[2];
sx q[2];
rz(2.810478) q[2];
rz(0.75774276) q[3];
sx q[3];
rz(-2.7527633) q[3];
sx q[3];
rz(-3.0537135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8452633) q[0];
sx q[0];
rz(-2.1305278) q[0];
sx q[0];
rz(2.9560126) q[0];
rz(-2.045385) q[1];
sx q[1];
rz(-2.9269693) q[1];
sx q[1];
rz(-1.4846444) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7916242) q[0];
sx q[0];
rz(-1.1317283) q[0];
sx q[0];
rz(-2.8429948) q[0];
x q[1];
rz(0.22557232) q[2];
sx q[2];
rz(-1.301287) q[2];
sx q[2];
rz(0.05664209) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2828196) q[1];
sx q[1];
rz(-0.38837896) q[1];
sx q[1];
rz(1.3057083) q[1];
rz(-pi) q[2];
rz(-0.0049403355) q[3];
sx q[3];
rz(-0.86588174) q[3];
sx q[3];
rz(-0.73965328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.93402702) q[2];
sx q[2];
rz(-0.88946122) q[2];
sx q[2];
rz(-0.55220848) q[2];
rz(2.3637555) q[3];
sx q[3];
rz(-0.85770291) q[3];
sx q[3];
rz(0.51789969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8778397) q[0];
sx q[0];
rz(-1.629566) q[0];
sx q[0];
rz(-1.4019479) q[0];
rz(-2.9539625) q[1];
sx q[1];
rz(-1.7788806) q[1];
sx q[1];
rz(2.3685041) q[1];
rz(0.4734584) q[2];
sx q[2];
rz(-0.31242328) q[2];
sx q[2];
rz(-1.7996126) q[2];
rz(-0.88541661) q[3];
sx q[3];
rz(-1.0049184) q[3];
sx q[3];
rz(0.64941209) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
