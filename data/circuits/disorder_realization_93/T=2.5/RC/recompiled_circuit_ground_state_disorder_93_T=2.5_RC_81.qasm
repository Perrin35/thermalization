OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.20004162) q[0];
sx q[0];
rz(3.3257421) q[0];
sx q[0];
rz(10.781547) q[0];
rz(1.4461203) q[1];
sx q[1];
rz(-1.5674011) q[1];
sx q[1];
rz(0.27624908) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8987792) q[0];
sx q[0];
rz(-1.587529) q[0];
sx q[0];
rz(-1.3748598) q[0];
x q[1];
rz(1.5803361) q[2];
sx q[2];
rz(-1.6838639) q[2];
sx q[2];
rz(2.557977) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3018926) q[1];
sx q[1];
rz(-3.0954315) q[1];
sx q[1];
rz(-0.34997527) q[1];
x q[2];
rz(2.802335) q[3];
sx q[3];
rz(-1.5087624) q[3];
sx q[3];
rz(1.3058939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9487298) q[2];
sx q[2];
rz(-2.7728045) q[2];
sx q[2];
rz(2.8249439) q[2];
rz(0.51490319) q[3];
sx q[3];
rz(-3.1334936) q[3];
sx q[3];
rz(-0.81075877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0136451) q[0];
sx q[0];
rz(-0.34334308) q[0];
sx q[0];
rz(-3.1105594) q[0];
rz(-1.7015999) q[1];
sx q[1];
rz(-0.80039918) q[1];
sx q[1];
rz(-1.7354785) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7714376) q[0];
sx q[0];
rz(-1.1509799) q[0];
sx q[0];
rz(0.19924723) q[0];
x q[1];
rz(1.3857394) q[2];
sx q[2];
rz(-0.20774797) q[2];
sx q[2];
rz(0.18321887) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1202218) q[1];
sx q[1];
rz(-1.7517266) q[1];
sx q[1];
rz(1.4814936) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7201028) q[3];
sx q[3];
rz(-0.79989767) q[3];
sx q[3];
rz(-0.61198583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7180353) q[2];
sx q[2];
rz(-2.9325298) q[2];
sx q[2];
rz(-2.5162856) q[2];
rz(1.5705671) q[3];
sx q[3];
rz(-2.8795241) q[3];
sx q[3];
rz(-0.03381332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5686947) q[0];
sx q[0];
rz(-0.53302065) q[0];
sx q[0];
rz(-1.1385338) q[0];
rz(-2.2472527) q[1];
sx q[1];
rz(-3.1411451) q[1];
sx q[1];
rz(2.0785887) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.759623) q[0];
sx q[0];
rz(-0.92734018) q[0];
sx q[0];
rz(-0.30610419) q[0];
x q[1];
rz(2.9350566) q[2];
sx q[2];
rz(-2.3828232) q[2];
sx q[2];
rz(1.1772924) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6394807) q[1];
sx q[1];
rz(-1.4993789) q[1];
sx q[1];
rz(1.4401431) q[1];
rz(2.9195905) q[3];
sx q[3];
rz(-1.3679805) q[3];
sx q[3];
rz(1.1403644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2669617) q[2];
sx q[2];
rz(-0.78445542) q[2];
sx q[2];
rz(2.0060914) q[2];
rz(1.9973495) q[3];
sx q[3];
rz(-0.091322986) q[3];
sx q[3];
rz(1.3967561) q[3];
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
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8046232) q[0];
sx q[0];
rz(-2.7380044) q[0];
sx q[0];
rz(2.023977) q[0];
rz(-3.0506253) q[1];
sx q[1];
rz(-0.03799835) q[1];
sx q[1];
rz(-1.631558) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9695797) q[0];
sx q[0];
rz(-1.0019127) q[0];
sx q[0];
rz(-2.8952231) q[0];
rz(2.6563866) q[2];
sx q[2];
rz(-2.910457) q[2];
sx q[2];
rz(2.427048) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9659057) q[1];
sx q[1];
rz(-1.693778) q[1];
sx q[1];
rz(1.2634274) q[1];
rz(-pi) q[2];
rz(-1.0311866) q[3];
sx q[3];
rz(-1.0301208) q[3];
sx q[3];
rz(0.074652076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9522004) q[2];
sx q[2];
rz(-0.11710937) q[2];
sx q[2];
rz(1.8899151) q[2];
rz(-0.39024726) q[3];
sx q[3];
rz(-3.1249505) q[3];
sx q[3];
rz(2.9515631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8613043) q[0];
sx q[0];
rz(-0.6328634) q[0];
sx q[0];
rz(1.2907775) q[0];
rz(3.1351008) q[1];
sx q[1];
rz(-0.26632729) q[1];
sx q[1];
rz(1.806462) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9296918) q[0];
sx q[0];
rz(-1.2294571) q[0];
sx q[0];
rz(1.7153738) q[0];
rz(-pi) q[1];
rz(-0.53369227) q[2];
sx q[2];
rz(-1.6385363) q[2];
sx q[2];
rz(2.2725032) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8346429) q[1];
sx q[1];
rz(-1.6942624) q[1];
sx q[1];
rz(2.929652) q[1];
x q[2];
rz(-1.174093) q[3];
sx q[3];
rz(-0.45214995) q[3];
sx q[3];
rz(-0.74739425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.73162762) q[2];
sx q[2];
rz(-0.16606398) q[2];
sx q[2];
rz(2.6483193) q[2];
rz(-1.0924245) q[3];
sx q[3];
rz(-3.137393) q[3];
sx q[3];
rz(-0.56001979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6888371) q[0];
sx q[0];
rz(-1.1863363) q[0];
sx q[0];
rz(1.3415927) q[0];
rz(-1.5635368) q[1];
sx q[1];
rz(-0.23461239) q[1];
sx q[1];
rz(2.2685307) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7811383) q[0];
sx q[0];
rz(-0.53780452) q[0];
sx q[0];
rz(-2.8685158) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6288343) q[2];
sx q[2];
rz(-1.1317694) q[2];
sx q[2];
rz(-0.88397938) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6775408) q[1];
sx q[1];
rz(-1.6455435) q[1];
sx q[1];
rz(1.4483553) q[1];
x q[2];
rz(-2.2544056) q[3];
sx q[3];
rz(-2.8018913) q[3];
sx q[3];
rz(1.8357133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0727498) q[2];
sx q[2];
rz(-1.9708956) q[2];
sx q[2];
rz(1.076131) q[2];
rz(2.7592646) q[3];
sx q[3];
rz(-0.040150661) q[3];
sx q[3];
rz(-0.7676777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1165987) q[0];
sx q[0];
rz(-3.0285663) q[0];
sx q[0];
rz(0.22079994) q[0];
rz(-0.91572541) q[1];
sx q[1];
rz(-2.974739) q[1];
sx q[1];
rz(-1.487026) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0512569) q[0];
sx q[0];
rz(-1.975734) q[0];
sx q[0];
rz(-1.9693518) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0887695) q[2];
sx q[2];
rz(-2.7332167) q[2];
sx q[2];
rz(2.0118304) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4055192) q[1];
sx q[1];
rz(-0.87784518) q[1];
sx q[1];
rz(0.939721) q[1];
rz(-0.12067699) q[3];
sx q[3];
rz(-1.7117705) q[3];
sx q[3];
rz(-1.3779061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23358341) q[2];
sx q[2];
rz(-0.073069409) q[2];
sx q[2];
rz(-2.0230587) q[2];
rz(-1.3633049) q[3];
sx q[3];
rz(-3.0647291) q[3];
sx q[3];
rz(-1.3479056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8175279) q[0];
sx q[0];
rz(-0.20453608) q[0];
sx q[0];
rz(1.2943785) q[0];
rz(1.4140465) q[1];
sx q[1];
rz(-3.0906446) q[1];
sx q[1];
rz(-2.9185413) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4572499) q[0];
sx q[0];
rz(-0.97829223) q[0];
sx q[0];
rz(2.1564583) q[0];
rz(-pi) q[1];
rz(1.4329918) q[2];
sx q[2];
rz(-1.3827822) q[2];
sx q[2];
rz(1.540198) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.053716252) q[1];
sx q[1];
rz(-1.1618335) q[1];
sx q[1];
rz(2.2323204) q[1];
rz(2.2659088) q[3];
sx q[3];
rz(-1.6611665) q[3];
sx q[3];
rz(2.4545057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.53905067) q[2];
sx q[2];
rz(-1.1844623) q[2];
sx q[2];
rz(0.43262532) q[2];
rz(-2.1116347) q[3];
sx q[3];
rz(-3.0410788) q[3];
sx q[3];
rz(2.4242134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3887727) q[0];
sx q[0];
rz(-2.6322375) q[0];
sx q[0];
rz(-0.21554047) q[0];
rz(1.132025) q[1];
sx q[1];
rz(-2.3339381) q[1];
sx q[1];
rz(-1.4707123) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12080321) q[0];
sx q[0];
rz(-1.258158) q[0];
sx q[0];
rz(2.1759869) q[0];
x q[1];
rz(-2.2756654) q[2];
sx q[2];
rz(-0.54958639) q[2];
sx q[2];
rz(1.7333502) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.52016693) q[1];
sx q[1];
rz(-1.5158719) q[1];
sx q[1];
rz(-1.5406233) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6588259) q[3];
sx q[3];
rz(-1.996187) q[3];
sx q[3];
rz(1.9874043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1620764) q[2];
sx q[2];
rz(-0.039204892) q[2];
sx q[2];
rz(-0.51852477) q[2];
rz(-2.3833852) q[3];
sx q[3];
rz(-0.77986217) q[3];
sx q[3];
rz(-1.6583748) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46490797) q[0];
sx q[0];
rz(-1.2228891) q[0];
sx q[0];
rz(1.2567047) q[0];
rz(-2.8769809) q[1];
sx q[1];
rz(-0.0022609641) q[1];
sx q[1];
rz(-1.1983926) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0459185) q[0];
sx q[0];
rz(-1.5642279) q[0];
sx q[0];
rz(1.571093) q[0];
rz(-pi) q[1];
rz(1.374515) q[2];
sx q[2];
rz(-0.33995865) q[2];
sx q[2];
rz(-0.46074552) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8598946) q[1];
sx q[1];
rz(-1.5571737) q[1];
sx q[1];
rz(2.9071684) q[1];
rz(-1.352574) q[3];
sx q[3];
rz(-0.17328158) q[3];
sx q[3];
rz(-0.97001433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4635224) q[2];
sx q[2];
rz(-0.0040201298) q[2];
sx q[2];
rz(3.1275952) q[2];
rz(2.8409581) q[3];
sx q[3];
rz(-0.0036573452) q[3];
sx q[3];
rz(-1.5716871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.32776253) q[0];
sx q[0];
rz(-1.7131085) q[0];
sx q[0];
rz(-1.499927) q[0];
rz(-0.035540237) q[1];
sx q[1];
rz(-2.7062369) q[1];
sx q[1];
rz(-2.8027262) q[1];
rz(0.69721194) q[2];
sx q[2];
rz(-1.7538191) q[2];
sx q[2];
rz(1.7303) q[2];
rz(-2.3382242) q[3];
sx q[3];
rz(-1.7163591) q[3];
sx q[3];
rz(1.8952445) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
