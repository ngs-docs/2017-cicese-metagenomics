===========================================================
Logging into your instance "in the cloud" (Windows version)
===========================================================

You'll need the .ppk file from the hackmd, and the hostname from the
spreadsheet, to log in.

Logging into your EC2 instance with Putty
=========================================

Open up putty, and enter your hostname into the Host Name box.

.. image:: images/win-putty-1.png
   :width: 50%

Now, go find the 'SSH' section and enter your ppk file (generated above
by puttygen).  Then select 'Open'.

.. image:: images/win-putty-2.png
   :width: 50%

Log in as "ubuntu".

.. image:: images/win-putty-3.png
   :width: 50%

Declare victory!

.. image:: images/win-putty-4.png
   :width: 50%

Here, you're logging in as user 'ubuntu' to the machine
'ec2-174-129-122-189.compute-1.amazonaws.com' using the authentication
key located in 'cicese.ppk' on your Desktop.

(Your key will be named 'cicese.ppk', and your hostname will be the
one in the spreadsheet.)

You should now see a text line that starts with something like
``ubuntu@ip-10-235-34-223:~$``.  You're in!  Now type::

   sudo bash
   cd /root

to switch into superuser mode (see: http://xkcd.com/149/) and go to your
home directory.

This is where the rest of the tutorials will start!

If you have Dropbox, you should now visit :doc:`installing-dropbox`.
To log out, type::

   exit
   logout

or just close the window.
