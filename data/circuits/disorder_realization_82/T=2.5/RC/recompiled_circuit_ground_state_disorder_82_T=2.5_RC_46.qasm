OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7656443) q[0];
sx q[0];
rz(-0.41416895) q[0];
sx q[0];
rz(2.3400657) q[0];
rz(-0.0030567788) q[1];
sx q[1];
rz(-0.86441511) q[1];
sx q[1];
rz(-0.094223082) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4163001) q[0];
sx q[0];
rz(-2.2803118) q[0];
sx q[0];
rz(0.50523357) q[0];
rz(-pi) q[1];
rz(1.8434486) q[2];
sx q[2];
rz(-1.6003719) q[2];
sx q[2];
rz(2.3675413) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4544782) q[1];
sx q[1];
rz(-1.2839437) q[1];
sx q[1];
rz(-2.0583365) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9889936) q[3];
sx q[3];
rz(-1.5729701) q[3];
sx q[3];
rz(-1.2945557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6039383) q[2];
sx q[2];
rz(-0.64017576) q[2];
sx q[2];
rz(-1.5113277) q[2];
rz(0.18621914) q[3];
sx q[3];
rz(-2.7112466) q[3];
sx q[3];
rz(-2.490624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6179825) q[0];
sx q[0];
rz(-1.1815434) q[0];
sx q[0];
rz(-0.13667662) q[0];
rz(-0.094820529) q[1];
sx q[1];
rz(-0.46351981) q[1];
sx q[1];
rz(-0.76900855) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.626132) q[0];
sx q[0];
rz(-0.40053408) q[0];
sx q[0];
rz(-0.35450165) q[0];
rz(0.42119512) q[2];
sx q[2];
rz(-0.73891089) q[2];
sx q[2];
rz(-0.74904672) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0598073) q[1];
sx q[1];
rz(-1.7331274) q[1];
sx q[1];
rz(-2.3947706) q[1];
rz(1.8212998) q[3];
sx q[3];
rz(-3.0022096) q[3];
sx q[3];
rz(-1.2047307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7424221) q[2];
sx q[2];
rz(-0.48226446) q[2];
sx q[2];
rz(1.4313618) q[2];
rz(-2.6509905) q[3];
sx q[3];
rz(-0.49121818) q[3];
sx q[3];
rz(-0.77559364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2630149) q[0];
sx q[0];
rz(-2.7624625) q[0];
sx q[0];
rz(-2.0468792) q[0];
rz(-1.8393983) q[1];
sx q[1];
rz(-0.98919386) q[1];
sx q[1];
rz(-3.1375695) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4152609) q[0];
sx q[0];
rz(-1.7021966) q[0];
sx q[0];
rz(0.8432998) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4039359) q[2];
sx q[2];
rz(-2.1263543) q[2];
sx q[2];
rz(-2.6839724) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0931172) q[1];
sx q[1];
rz(-2.8826536) q[1];
sx q[1];
rz(0.11872752) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4278042) q[3];
sx q[3];
rz(-1.3828438) q[3];
sx q[3];
rz(-0.97602188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5502988) q[2];
sx q[2];
rz(-0.58612263) q[2];
sx q[2];
rz(-2.6934521) q[2];
rz(0.48629931) q[3];
sx q[3];
rz(-0.86217642) q[3];
sx q[3];
rz(-1.1264616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9124741) q[0];
sx q[0];
rz(-1.7121226) q[0];
sx q[0];
rz(-2.4718156) q[0];
rz(-1.229333) q[1];
sx q[1];
rz(-0.45893097) q[1];
sx q[1];
rz(0.545048) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16825039) q[0];
sx q[0];
rz(-1.1118044) q[0];
sx q[0];
rz(-2.5516872) q[0];
rz(-0.82370211) q[2];
sx q[2];
rz(-0.81497008) q[2];
sx q[2];
rz(-2.4832249) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.5601756) q[1];
sx q[1];
rz(-0.22318527) q[1];
sx q[1];
rz(1.8471922) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0470951) q[3];
sx q[3];
rz(-2.3102488) q[3];
sx q[3];
rz(0.0090713105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5547319) q[2];
sx q[2];
rz(-1.6394337) q[2];
sx q[2];
rz(0.16793212) q[2];
rz(1.933291) q[3];
sx q[3];
rz(-0.30253634) q[3];
sx q[3];
rz(-1.6944073) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5912882) q[0];
sx q[0];
rz(-1.8809603) q[0];
sx q[0];
rz(3.1377129) q[0];
rz(0.30715352) q[1];
sx q[1];
rz(-0.75633621) q[1];
sx q[1];
rz(-0.53877962) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2870432) q[0];
sx q[0];
rz(-1.3417305) q[0];
sx q[0];
rz(-1.4255217) q[0];
rz(-pi) q[1];
rz(1.4319375) q[2];
sx q[2];
rz(-1.9966239) q[2];
sx q[2];
rz(2.6151997) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6314124) q[1];
sx q[1];
rz(-2.4459527) q[1];
sx q[1];
rz(2.5744777) q[1];
rz(-pi) q[2];
rz(-0.73154463) q[3];
sx q[3];
rz(-1.9172102) q[3];
sx q[3];
rz(-2.9820095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4469845) q[2];
sx q[2];
rz(-1.0993404) q[2];
sx q[2];
rz(1.851932) q[2];
rz(1.5223711) q[3];
sx q[3];
rz(-2.3963942) q[3];
sx q[3];
rz(1.4815909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.30037844) q[0];
sx q[0];
rz(-2.2024246) q[0];
sx q[0];
rz(0.10264957) q[0];
rz(2.4094021) q[1];
sx q[1];
rz(-0.79473549) q[1];
sx q[1];
rz(-0.57714677) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6711376) q[0];
sx q[0];
rz(-1.5671434) q[0];
sx q[0];
rz(1.5421583) q[0];
x q[1];
rz(-2.5546165) q[2];
sx q[2];
rz(-1.7086626) q[2];
sx q[2];
rz(-1.1816813) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0476372) q[1];
sx q[1];
rz(-3.006662) q[1];
sx q[1];
rz(0.62420242) q[1];
rz(-pi) q[2];
rz(-1.243337) q[3];
sx q[3];
rz(-2.2994882) q[3];
sx q[3];
rz(3.1407243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.010043667) q[2];
sx q[2];
rz(-0.71655822) q[2];
sx q[2];
rz(1.6183759) q[2];
rz(-3.0736382) q[3];
sx q[3];
rz(-0.32502919) q[3];
sx q[3];
rz(2.3006191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1286569) q[0];
sx q[0];
rz(-2.1067297) q[0];
sx q[0];
rz(-0.77142429) q[0];
rz(2.6029288) q[1];
sx q[1];
rz(-1.795105) q[1];
sx q[1];
rz(1.0661941) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26771256) q[0];
sx q[0];
rz(-1.6475878) q[0];
sx q[0];
rz(1.3174873) q[0];
rz(-pi) q[1];
rz(1.0577002) q[2];
sx q[2];
rz(-1.7745681) q[2];
sx q[2];
rz(0.17901006) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.297841) q[1];
sx q[1];
rz(-0.21756309) q[1];
sx q[1];
rz(3.080653) q[1];
x q[2];
rz(2.4769884) q[3];
sx q[3];
rz(-0.41575894) q[3];
sx q[3];
rz(-0.98008728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44925877) q[2];
sx q[2];
rz(-0.80654311) q[2];
sx q[2];
rz(2.6991357) q[2];
rz(-0.19065204) q[3];
sx q[3];
rz(-2.4176044) q[3];
sx q[3];
rz(2.382522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8964748) q[0];
sx q[0];
rz(-2.9560095) q[0];
sx q[0];
rz(1.8316487) q[0];
rz(-2.7438121) q[1];
sx q[1];
rz(-2.3186627) q[1];
sx q[1];
rz(-1.3936874) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8373972) q[0];
sx q[0];
rz(-2.1050354) q[0];
sx q[0];
rz(-0.5128575) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2692503) q[2];
sx q[2];
rz(-1.796495) q[2];
sx q[2];
rz(1.7762453) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2893148) q[1];
sx q[1];
rz(-0.56486579) q[1];
sx q[1];
rz(-2.5603065) q[1];
rz(-pi) q[2];
rz(1.0839127) q[3];
sx q[3];
rz(-2.9176035) q[3];
sx q[3];
rz(-1.4518224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1132249) q[2];
sx q[2];
rz(-0.71030474) q[2];
sx q[2];
rz(-0.073471546) q[2];
rz(2.4469589) q[3];
sx q[3];
rz(-1.8690542) q[3];
sx q[3];
rz(-1.0276444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46591127) q[0];
sx q[0];
rz(-2.9407192) q[0];
sx q[0];
rz(-2.1821816) q[0];
rz(-0.77964669) q[1];
sx q[1];
rz(-1.0028853) q[1];
sx q[1];
rz(0.27228212) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5648266) q[0];
sx q[0];
rz(-1.5656002) q[0];
sx q[0];
rz(-1.5840498) q[0];
rz(-pi) q[1];
rz(-1.9302888) q[2];
sx q[2];
rz(-0.73795107) q[2];
sx q[2];
rz(0.28970813) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.694937) q[1];
sx q[1];
rz(-2.4866101) q[1];
sx q[1];
rz(-2.2769364) q[1];
rz(-pi) q[2];
rz(-2.394478) q[3];
sx q[3];
rz(-2.0339587) q[3];
sx q[3];
rz(-0.85827352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.73725629) q[2];
sx q[2];
rz(-1.5055089) q[2];
sx q[2];
rz(0.20757248) q[2];
rz(3.0370144) q[3];
sx q[3];
rz(-0.50555491) q[3];
sx q[3];
rz(1.6149717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.680147) q[0];
sx q[0];
rz(-1.2274281) q[0];
sx q[0];
rz(1.0567868) q[0];
rz(2.5959004) q[1];
sx q[1];
rz(-0.97815198) q[1];
sx q[1];
rz(-0.32870764) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0595039) q[0];
sx q[0];
rz(-0.57841136) q[0];
sx q[0];
rz(0.93579458) q[0];
rz(0.64735554) q[2];
sx q[2];
rz(-0.090687625) q[2];
sx q[2];
rz(-2.2272404) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.94241566) q[1];
sx q[1];
rz(-0.34245447) q[1];
sx q[1];
rz(-1.5210995) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28139754) q[3];
sx q[3];
rz(-0.79165047) q[3];
sx q[3];
rz(0.36330128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9823965) q[2];
sx q[2];
rz(-3.0089162) q[2];
sx q[2];
rz(2.0170434) q[2];
rz(-3.0647965) q[3];
sx q[3];
rz(-2.7087637) q[3];
sx q[3];
rz(2.2176149) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15976739) q[0];
sx q[0];
rz(-1.2564909) q[0];
sx q[0];
rz(1.5782574) q[0];
rz(-2.3253597) q[1];
sx q[1];
rz(-1.4030133) q[1];
sx q[1];
rz(-1.0907008) q[1];
rz(2.6098567) q[2];
sx q[2];
rz(-0.81960631) q[2];
sx q[2];
rz(-1.3255957) q[2];
rz(2.7362551) q[3];
sx q[3];
rz(-1.574285) q[3];
sx q[3];
rz(-2.3226572) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
