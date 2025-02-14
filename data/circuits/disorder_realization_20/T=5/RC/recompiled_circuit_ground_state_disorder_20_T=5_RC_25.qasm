OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5372758) q[0];
sx q[0];
rz(-0.24157) q[0];
sx q[0];
rz(-2.8085652) q[0];
rz(2.0060519) q[1];
sx q[1];
rz(-0.82692868) q[1];
sx q[1];
rz(2.4976322) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4304023) q[0];
sx q[0];
rz(-1.2997753) q[0];
sx q[0];
rz(-2.2008116) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7276554) q[2];
sx q[2];
rz(-1.2741977) q[2];
sx q[2];
rz(3.0576884) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8713069) q[1];
sx q[1];
rz(-1.3340923) q[1];
sx q[1];
rz(-0.38856296) q[1];
x q[2];
rz(-1.1584362) q[3];
sx q[3];
rz(-1.8915911) q[3];
sx q[3];
rz(-3.0123346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.48646271) q[2];
sx q[2];
rz(-2.6772406) q[2];
sx q[2];
rz(-2.4008524) q[2];
rz(1.5517392) q[3];
sx q[3];
rz(-2.430075) q[3];
sx q[3];
rz(2.5488502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43316677) q[0];
sx q[0];
rz(-1.1335224) q[0];
sx q[0];
rz(-3.0246227) q[0];
rz(2.6372657) q[1];
sx q[1];
rz(-1.802899) q[1];
sx q[1];
rz(0.28409827) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5829651) q[0];
sx q[0];
rz(-0.57218781) q[0];
sx q[0];
rz(-1.7250604) q[0];
x q[1];
rz(-2.5710377) q[2];
sx q[2];
rz(-2.2584256) q[2];
sx q[2];
rz(-0.90435435) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.85728474) q[1];
sx q[1];
rz(-1.6928288) q[1];
sx q[1];
rz(-0.034842592) q[1];
x q[2];
rz(-2.8900364) q[3];
sx q[3];
rz(-1.8402035) q[3];
sx q[3];
rz(-0.081307383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8154907) q[2];
sx q[2];
rz(-0.88560605) q[2];
sx q[2];
rz(-0.76078129) q[2];
rz(-0.86959362) q[3];
sx q[3];
rz(-1.2660916) q[3];
sx q[3];
rz(0.45309711) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3608383) q[0];
sx q[0];
rz(-2.0212845) q[0];
sx q[0];
rz(-0.25217062) q[0];
rz(-1.699327) q[1];
sx q[1];
rz(-2.1863054) q[1];
sx q[1];
rz(2.7853277) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7423681) q[0];
sx q[0];
rz(-0.059905298) q[0];
sx q[0];
rz(-1.674288) q[0];
rz(-pi) q[1];
rz(2.5926931) q[2];
sx q[2];
rz(-1.2147012) q[2];
sx q[2];
rz(-2.9402972) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2846996) q[1];
sx q[1];
rz(-2.2027122) q[1];
sx q[1];
rz(-0.17669682) q[1];
rz(-pi) q[2];
rz(-1.3144794) q[3];
sx q[3];
rz(-1.5689092) q[3];
sx q[3];
rz(-2.6607115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0226655) q[2];
sx q[2];
rz(-1.3902731) q[2];
sx q[2];
rz(1.5125795) q[2];
rz(3.1318943) q[3];
sx q[3];
rz(-1.9111948) q[3];
sx q[3];
rz(-2.5201918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-1.6464722) q[0];
sx q[0];
rz(-0.54238129) q[0];
sx q[0];
rz(0.25076184) q[0];
rz(1.3465025) q[1];
sx q[1];
rz(-0.52918068) q[1];
sx q[1];
rz(1.4432602) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0899857) q[0];
sx q[0];
rz(-2.3653226) q[0];
sx q[0];
rz(-1.8818186) q[0];
x q[1];
rz(-2.2002831) q[2];
sx q[2];
rz(-2.513859) q[2];
sx q[2];
rz(2.8986487) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52407284) q[1];
sx q[1];
rz(-1.0403324) q[1];
sx q[1];
rz(2.4071715) q[1];
rz(-pi) q[2];
x q[2];
rz(0.26224995) q[3];
sx q[3];
rz(-1.831154) q[3];
sx q[3];
rz(0.32338705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.5556339) q[2];
sx q[2];
rz(-1.1052174) q[2];
sx q[2];
rz(2.8982758) q[2];
rz(-2.5569052) q[3];
sx q[3];
rz(-0.57716113) q[3];
sx q[3];
rz(-3.1134591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1481767) q[0];
sx q[0];
rz(-2.7758444) q[0];
sx q[0];
rz(-0.17955968) q[0];
rz(-1.2252294) q[1];
sx q[1];
rz(-1.4045249) q[1];
sx q[1];
rz(-2.6833351) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1711463) q[0];
sx q[0];
rz(-1.2752646) q[0];
sx q[0];
rz(-1.4895205) q[0];
x q[1];
rz(1.304428) q[2];
sx q[2];
rz(-0.86037105) q[2];
sx q[2];
rz(-1.2183587) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.51650713) q[1];
sx q[1];
rz(-1.436215) q[1];
sx q[1];
rz(-0.075304042) q[1];
x q[2];
rz(1.1510994) q[3];
sx q[3];
rz(-2.3956798) q[3];
sx q[3];
rz(1.0585566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3852343) q[2];
sx q[2];
rz(-0.42432722) q[2];
sx q[2];
rz(2.2198086) q[2];
rz(-1.8900169) q[3];
sx q[3];
rz(-1.7332964) q[3];
sx q[3];
rz(-3.0799227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.09403041) q[0];
sx q[0];
rz(-0.91218364) q[0];
sx q[0];
rz(-0.14866522) q[0];
rz(0.31691638) q[1];
sx q[1];
rz(-2.6628559) q[1];
sx q[1];
rz(-0.61683488) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4267093) q[0];
sx q[0];
rz(-1.6322109) q[0];
sx q[0];
rz(-1.6493357) q[0];
x q[1];
rz(-2.0111175) q[2];
sx q[2];
rz(-1.2424349) q[2];
sx q[2];
rz(0.87000123) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9636391) q[1];
sx q[1];
rz(-2.511189) q[1];
sx q[1];
rz(-2.0042581) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1149949) q[3];
sx q[3];
rz(-0.57145703) q[3];
sx q[3];
rz(-2.8570017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2769015) q[2];
sx q[2];
rz(-1.428705) q[2];
sx q[2];
rz(-3.1332664) q[2];
rz(-0.62638038) q[3];
sx q[3];
rz(-2.4255987) q[3];
sx q[3];
rz(-0.00055073784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.83547) q[0];
sx q[0];
rz(-1.1518814) q[0];
sx q[0];
rz(-0.48208153) q[0];
rz(-3.112402) q[1];
sx q[1];
rz(-1.2960641) q[1];
sx q[1];
rz(0.76404244) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8523435) q[0];
sx q[0];
rz(-1.2004392) q[0];
sx q[0];
rz(-0.072288805) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5381728) q[2];
sx q[2];
rz(-2.4815686) q[2];
sx q[2];
rz(0.78885022) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5611539) q[1];
sx q[1];
rz(-2.517546) q[1];
sx q[1];
rz(1.2777722) q[1];
rz(-pi) q[2];
rz(0.0079844012) q[3];
sx q[3];
rz(-1.1698876) q[3];
sx q[3];
rz(2.8454091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.46334106) q[2];
sx q[2];
rz(-1.8222787) q[2];
sx q[2];
rz(-0.48416644) q[2];
rz(0.68228996) q[3];
sx q[3];
rz(-2.4874918) q[3];
sx q[3];
rz(-3.0068523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.057137) q[0];
sx q[0];
rz(-2.3637922) q[0];
sx q[0];
rz(-1.2917668) q[0];
rz(-2.6765587) q[1];
sx q[1];
rz(-0.51883042) q[1];
sx q[1];
rz(0.30716392) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43873337) q[0];
sx q[0];
rz(-0.68638681) q[0];
sx q[0];
rz(0.86875654) q[0];
rz(-pi) q[1];
rz(-2.6382006) q[2];
sx q[2];
rz(-0.97676986) q[2];
sx q[2];
rz(-0.6436178) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0867566) q[1];
sx q[1];
rz(-1.0030748) q[1];
sx q[1];
rz(-2.9556429) q[1];
rz(2.516496) q[3];
sx q[3];
rz(-1.5872883) q[3];
sx q[3];
rz(1.0436637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.18053599) q[2];
sx q[2];
rz(-0.44027105) q[2];
sx q[2];
rz(2.4214936) q[2];
rz(-0.85047203) q[3];
sx q[3];
rz(-1.1668147) q[3];
sx q[3];
rz(-1.7299962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9367323) q[0];
sx q[0];
rz(-0.16848773) q[0];
sx q[0];
rz(2.4334461) q[0];
rz(1.6234966) q[1];
sx q[1];
rz(-0.93815362) q[1];
sx q[1];
rz(2.785397) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4771381) q[0];
sx q[0];
rz(-1.150048) q[0];
sx q[0];
rz(-0.53945213) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5237066) q[2];
sx q[2];
rz(-1.5466585) q[2];
sx q[2];
rz(1.2302421) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1005786) q[1];
sx q[1];
rz(-1.6008401) q[1];
sx q[1];
rz(-2.492659) q[1];
rz(-pi) q[2];
rz(1.1061927) q[3];
sx q[3];
rz(-1.2038017) q[3];
sx q[3];
rz(1.3657686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0922962) q[2];
sx q[2];
rz(-2.3342817) q[2];
sx q[2];
rz(0.88579196) q[2];
rz(0.93585912) q[3];
sx q[3];
rz(-1.1805308) q[3];
sx q[3];
rz(0.22127557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43779272) q[0];
sx q[0];
rz(-3.0942823) q[0];
sx q[0];
rz(1.4495151) q[0];
rz(-2.119078) q[1];
sx q[1];
rz(-2.6578564) q[1];
sx q[1];
rz(-1.6922916) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.463844) q[0];
sx q[0];
rz(-0.41526702) q[0];
sx q[0];
rz(1.4474611) q[0];
rz(-pi) q[1];
rz(-0.42745356) q[2];
sx q[2];
rz(-1.1910254) q[2];
sx q[2];
rz(-2.3878271) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.199618) q[1];
sx q[1];
rz(-1.7829547) q[1];
sx q[1];
rz(2.1320264) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14808296) q[3];
sx q[3];
rz(-1.8335153) q[3];
sx q[3];
rz(-2.1488291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8861822) q[2];
sx q[2];
rz(-1.1899199) q[2];
sx q[2];
rz(0.6748684) q[2];
rz(-2.1700962) q[3];
sx q[3];
rz(-2.4392305) q[3];
sx q[3];
rz(1.491588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7991199) q[0];
sx q[0];
rz(-1.5958888) q[0];
sx q[0];
rz(2.2842443) q[0];
rz(-2.9091861) q[1];
sx q[1];
rz(-2.1288165) q[1];
sx q[1];
rz(-1.7850599) q[1];
rz(-2.1017553) q[2];
sx q[2];
rz(-1.2086443) q[2];
sx q[2];
rz(2.3966387) q[2];
rz(2.2994249) q[3];
sx q[3];
rz(-1.5123868) q[3];
sx q[3];
rz(2.2342891) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
