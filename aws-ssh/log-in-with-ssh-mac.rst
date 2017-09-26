=======================================================
Logging into your instance "in the cloud" (Mac version)
=======================================================

OK, so we've given you a running computer. How do you access it?

The two thing you'll need are the network name of your computer, which
you can get from the spreadsheet linked to in the hackmd; and the
'cicese.pem' file that we've given you on the hackmd.

Copy the name, and connect to that computer with ssh under the username
'ubuntu', as follows:

First, move your private key file from wherever you downloaded it onto
your Desktop.

Next, start Terminal (in Applications... Utilities...) and type::

  chmod og-rwx ~/Desktop/cicese.pem

to set the permissions on the private key file to "closed to others".

Then type::

  ssh -i ~/Desktop/cicese.pem ubuntu@ec2-???-???-???-???.compute-1.amazonaws.com

(but you have to replace the stuff after the '@' sign with the name of the host).

Here, you're logging in as user 'ubuntu' to the machine
'ec2-174-129-122-189.compute-1.amazonaws.com' using the authentication
key located in 'cicese.pem' on your Desktop.

You should now see a text line that starts with something like
``ubuntu@ip-10-235-34-223:~$``.  You're in!  Now type::

   sudo bash
   cd /root

to switch into superuser mode (see: http://xkcd.com/149/) and go to your
home directory.

This is where the rest of the tutorials will start!

To log out, type::

   exit
   logout

or just close the window.
